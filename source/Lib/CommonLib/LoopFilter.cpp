/* -----------------------------------------------------------------------------
The copyright in this software is being made available under the BSD
License, included below. No patent rights, trademark rights and/or
other Intellectual Property Rights other than the copyrights concerning
the Software are granted under this license.

For any license concerning other Intellectual Property rights than the software,
especially patent licenses, a separate Agreement needs to be closed.
For more information please contact:

Fraunhofer Heinrich Hertz Institute
Einsteinufer 37
10587 Berlin, Germany
www.hhi.fraunhofer.de/vvc
vvc@hhi.fraunhofer.de

Copyright (c) 2018-2020, Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Fraunhofer nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.


------------------------------------------------------------------------------------------- */

/** \file     LoopFilter.cpp
    \brief    deblocking filter
*/

#include "LoopFilter.h"
#include "Slice.h"
#include "Mv.h"
#include "Unit.h"
#include "UnitTools.h"
#include "UnitPartitioner.h"
#include "dtrace_codingstruct.h"
#include "dtrace_buffer.h"
#include "CommonLib/TimeProfiler.h"

#include "Quant.h"

#ifdef TARGET_SIMD_X86
#  include "CommonDefX86.h"
#endif

//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Constants
// ====================================================================================================================

#define DEBLOCK_SMALLEST_BLOCK 8

#define DEFAULT_INTRA_TC_OFFSET 2  ///< Default intra TC offset

// ====================================================================================================================
// Tables
// ====================================================================================================================

const uint16_t LoopFilter::sm_tcTable[MAX_QP + 1 + DEFAULT_INTRA_TC_OFFSET] = {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   0,   0,   0,   3,   4,   4,   4,
    4,  5,  5,  5,  5,  7,  7,  8,  9,  10,  10,  11,  13,  14,  15,  17,  19,  21,  24,  25,  29,  33,
    36, 41, 45, 51, 57, 64, 71, 80, 89, 100, 112, 125, 141, 157, 177, 198, 222, 250, 280, 314, 352, 395};
const uint8_t LoopFilter::sm_betaTable[MAX_QP + 1] = {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                                                      6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 24,
                                                      26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56,
                                                      58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88};

// ====================================================================================================================
// utility functions
// ====================================================================================================================

template <bool isChromaHorCTBBoundary, typename T>
static inline int xCalcDP(const T* piSrc, const ptrdiff_t iOffset) {
  if (isChromaHorCTBBoundary) {
    return abs(piSrc[-iOffset * 2] - 2 * piSrc[-iOffset * 2] + piSrc[-iOffset]);
  } else {
    return abs(piSrc[-iOffset * 3] - 2 * piSrc[-iOffset * 2] + piSrc[-iOffset]);
  }
}

template <typename T>
static inline int xCalcDQ(const T* piSrc, const ptrdiff_t iOffset) {
  return abs(piSrc[0] - 2 * piSrc[iOffset] + piSrc[iOffset * 2]);
}

template <typename T>
static inline bool xUseStrongFiltering(const T* piSrc, const ptrdiff_t iOffset, const int d, const int beta,
                                       const int tc, bool sidePisLarge = false, bool sideQisLarge = false,
                                       int maxFilterLengthP = 7, int maxFilterLengthQ = 7,
                                       bool isChromaHorCTBBoundary = false) {
  const T m3 = piSrc[-1 * iOffset];
  const T m4 = piSrc[0];

  if (!(d < (beta >> 2) && abs(m3 - m4) < ((tc * 5 + 1) >> 1))) return false;

  const T m0 = piSrc[-4 * iOffset];
  const T m7 = piSrc[3 * iOffset];

  const T m2 = piSrc[-iOffset * 2];
  int sp3 = abs(m0 - m3);
  if (isChromaHorCTBBoundary) {
    sp3 = abs(m2 - m3);
  }
  int sq3 = abs(m7 - m4);
  const int d_strong = sp3 + sq3;

  if (sidePisLarge || sideQisLarge) {
    if (sidePisLarge) {
      const T mP4 = piSrc[-iOffset * maxFilterLengthP - iOffset];
      if (maxFilterLengthP == 7) {
        const T mP5 = piSrc[-iOffset * 5];
        const T mP6 = piSrc[-iOffset * 6];
        const T mP7 = piSrc[-iOffset * 7];
        sp3 = sp3 + abs(mP5 - mP6 - mP7 + mP4);
      }
      sp3 = (sp3 + abs(m0 - mP4) + 1) >> 1;
    }
    if (sideQisLarge) {
      const T m11 = piSrc[iOffset * maxFilterLengthQ];
      if (maxFilterLengthQ == 7) {
        const T m8 = piSrc[iOffset * 4];
        const T m9 = piSrc[iOffset * 5];
        const T m10 = piSrc[iOffset * 6];
        sq3 = sq3 + abs(m8 - m9 - m10 + m11);
      }
      sq3 = (sq3 + abs(m11 - m7) + 1) >> 1;
    }
    return ((sp3 + sq3) < (beta * 3 >> 5)) &&
           (d < (beta >> 4));
  } else {
    return d_strong < (beta >> 3);
  }
}

static bool isAvailable(const CodingUnit& cu, const CodingUnit& cu2, const bool bEnforceSliceRestriction,
                        const bool bEnforceTileRestriction, const bool bEnforceSubPicRestriction) {
  return ((!bEnforceSliceRestriction || CU::isSameSlice(cu, cu2)) &&
          (!bEnforceTileRestriction || CU::isSameTile(cu, cu2)) &&
          (!bEnforceSubPicRestriction || CU::isSameSubPic(cu, cu2)));
}

static int xCalBsSameRefCore(const MotionInfo& miP, const MotionInfo& miQ, const Picture* piRefP0,
                             const Picture* piRefP1, const Picture* piRefQ0, const Picture* piRefQ1) {
  static constexpr int nThreshold = (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 1;
  int uiBs;

  Mv mvP[2] = {{0, 0}, {0, 0}}, mvQ[2] = {{0, 0}, {0, 0}};

  if (0 <= miP.refIdx[0]) {
    mvP[0] = miP.mv[0];
  }
  if (0 <= miP.refIdx[1]) {
    mvP[1] = miP.mv[1];
  }
  if (0 <= miQ.refIdx[0]) {
    mvQ[0] = miQ.mv[0];
  }
  if (0 <= miQ.refIdx[1]) {
    mvQ[1] = miQ.mv[1];
  }
  if (piRefP0 != piRefP1) {  // Different L0 & L1
    if (piRefP0 == piRefQ0) {
      uiBs = ((abs(mvQ[0].getHor() - mvP[0].getHor()) >= nThreshold) ||
              (abs(mvQ[0].getVer() - mvP[0].getVer()) >= nThreshold) ||
              (abs(mvQ[1].getHor() - mvP[1].getHor()) >= nThreshold) ||
              (abs(mvQ[1].getVer() - mvP[1].getVer()) >= nThreshold))
                 ? 1
                 : 0;
    } else {
      uiBs = ((abs(mvQ[1].getHor() - mvP[0].getHor()) >= nThreshold) ||
              (abs(mvQ[1].getVer() - mvP[0].getVer()) >= nThreshold) ||
              (abs(mvQ[0].getHor() - mvP[1].getHor()) >= nThreshold) ||
              (abs(mvQ[0].getVer() - mvP[1].getVer()) >= nThreshold))
                 ? 1
                 : 0;
    }
  } else {  // Same L0 & L1
    uiBs = ((abs(mvQ[0].getHor() - mvP[0].getHor()) >= nThreshold) ||
            (abs(mvQ[0].getVer() - mvP[0].getVer()) >= nThreshold) ||
            (abs(mvQ[1].getHor() - mvP[1].getHor()) >= nThreshold) ||
            (abs(mvQ[1].getVer() - mvP[1].getVer()) >= nThreshold)) &&
                   ((abs(mvQ[1].getHor() - mvP[0].getHor()) >= nThreshold) ||
                    (abs(mvQ[1].getVer() - mvP[0].getVer()) >= nThreshold) ||
                    (abs(mvQ[0].getHor() - mvP[1].getHor()) >= nThreshold) ||
                    (abs(mvQ[0].getVer() - mvP[1].getVer()) >= nThreshold))
               ? 1
               : 0;
  }

  return uiBs;
}

#define BsSet(val, compIdx) ((val) << ((compIdx) << 1))
#define BsGet(val, compIdx) (((val) >> ((compIdx) << 1)) & 3)

static const int dbCoeffs7[7] = {59, 50, 41, 32, 23, 14, 5};
static const int dbCoeffs5[5] = {58, 45, 32, 19, 6};
static const int dbCoeffs3[3] = {53, 32, 11};

template <typename TSrcDst>
static inline void xBilinearFilter(TSrcDst* srcP, TSrcDst* srcQ, ptrdiff_t offset, int refMiddle, int refP, int refQ,
                                   int numberPSide, int numberQSide, const int* dbCoeffsP, const int* dbCoeffsQ,
                                   int tc) {
  int src;
  const char tc7[7] = {6, 5, 4, 3, 2, 1, 1};
  const char tc3[3] = {6, 4, 2};
  const char* tcP = (numberPSide == 3) ? tc3 : tc7;
  const char* tcQ = (numberQSide == 3) ? tc3 : tc7;

  for (int pos = 0; pos < numberPSide; pos++) {
    src = srcP[-offset * pos];
    int cvalue = (tc * tcP[pos]) >> 1;
    srcP[-offset * pos] =
        Clip3(src - cvalue, src + cvalue, ((refMiddle * dbCoeffsP[pos] + refP * (64 - dbCoeffsP[pos]) + 32) >> 6));
  }

  for (int pos = 0; pos < numberQSide; pos++) {
    src = srcQ[offset * pos];
    int cvalue = (tc * tcQ[pos]) >> 1;
    srcQ[offset * pos] =
        Clip3(src - cvalue, src + cvalue, ((refMiddle * dbCoeffsQ[pos] + refQ * (64 - dbCoeffsQ[pos]) + 32) >> 6));
  }
}

template <typename TSrcDst>
void xFilteringPandQCore(Pel* _src, ptrdiff_t step, const ptrdiff_t offset, int numberPSide, int numberQSide, int tc) {
  TSrcDst* src = reinterpret_cast<TSrcDst*>(_src);
  CHECK(numberPSide <= 3 && numberQSide <= 3, "Short filtering in long filtering function");
  CHECK(numberPSide != 3 && numberPSide != 5 && numberPSide != 7, "invalid numberPSide");
  CHECK(numberQSide != 3 && numberQSide != 5 && numberQSide != 7, "invalid numberQSide");

  const int* dbCoeffsP = numberPSide == 7 ? dbCoeffs7 : (numberPSide == 5) ? dbCoeffs5 : dbCoeffs3;
  const int* dbCoeffsQ = numberQSide == 7 ? dbCoeffs7 : (numberQSide == 5) ? dbCoeffs5 : dbCoeffs3;

  for (int i = 0; i < DEBLOCK_SMALLEST_BLOCK / 2; i++) {
    TSrcDst* srcP = src + step * i - offset;
    TSrcDst* srcQ = src + step * i;

    int refP = 0;
    int refQ = 0;
    int refMiddle = 0;

    switch (numberPSide) {
      case 7:
        refP = (srcP[-6 * offset] + srcP[-7 * offset] + 1) >> 1;
        break;
      case 3:
        refP = (srcP[-2 * offset] + srcP[-3 * offset] + 1) >> 1;
        break;
      case 5:
        refP = (srcP[-4 * offset] + srcP[-5 * offset] + 1) >> 1;
        break;
    }

    switch (numberQSide) {
      case 7:
        refQ = (srcQ[6 * offset] + srcQ[7 * offset] + 1) >> 1;
        break;
      case 3:
        refQ = (srcQ[2 * offset] + srcQ[3 * offset] + 1) >> 1;
        break;
      case 5:
        refQ = (srcQ[4 * offset] + srcQ[5 * offset] + 1) >> 1;
        break;
    }

    if (numberPSide == numberQSide) {
      if (numberPSide == 5) {
        refMiddle = (2 * (srcP[0] + srcQ[0] + srcP[-offset] + srcQ[offset] + srcP[-2 * offset] + srcQ[2 * offset]) +
                     srcP[-3 * offset] + srcQ[3 * offset] + srcP[-4 * offset] + srcQ[4 * offset] + 8) >>
                    4;
      } else {
        refMiddle = (2 * (srcP[0] + srcQ[0]) + srcP[-offset] + srcQ[offset] + srcP[-2 * offset] + srcQ[2 * offset] +
                     srcP[-3 * offset] + srcQ[3 * offset] + srcP[-4 * offset] + srcQ[4 * offset] + srcP[-5 * offset] +
                     srcQ[5 * offset] + +srcP[-6 * offset] + srcQ[6 * offset] + 8) >>
                    4;
      }
    } else {
      TSrcDst* srcPt = srcP;
      TSrcDst* srcQt = srcQ;
      ptrdiff_t offsetP = -offset;
      ptrdiff_t offsetQ = offset;

      int newNumberQSide = numberQSide;
      int newNumberPSide = numberPSide;

      if (numberQSide > numberPSide) {
        std::swap(srcPt, srcQt);
        std::swap(offsetP, offsetQ);
        newNumberQSide = numberPSide;
        newNumberPSide = numberQSide;
      }

      if (newNumberPSide == 7 && newNumberQSide == 5) {
        refMiddle = (2 * (srcP[0] + srcQ[0] + srcP[-offset] + srcQ[offset]) + srcP[-2 * offset] + srcQ[2 * offset] +
                     srcP[-3 * offset] + srcQ[3 * offset] + srcP[-4 * offset] + srcQ[4 * offset] + srcP[-5 * offset] +
                     srcQ[5 * offset] + 8) >>
                    4;
      } else if (newNumberPSide == 7 && newNumberQSide == 3) {
        refMiddle = (2 * (srcPt[0] + srcQt[0]) + srcQt[0] + 2 * (srcQt[offsetQ] + srcQt[2 * offsetQ]) + srcPt[offsetP] +
                     srcQt[offsetQ] + srcPt[2 * offsetP] + srcPt[3 * offsetP] + srcPt[4 * offsetP] +
                     srcPt[5 * offsetP] + srcPt[6 * offsetP] + 8) >>
                    4;
      } else {  // if (newNumberPSide == 5 && newNumberQSide == 3)
        refMiddle = (srcP[0] + srcQ[0] + srcP[-offset] + srcQ[offset] + srcP[-2 * offset] + srcQ[2 * offset] +
                     srcP[-3 * offset] + srcQ[3 * offset] + 4) >>
                    3;
      }
    }

    xBilinearFilter<TSrcDst>(srcP, srcQ, offset, refMiddle, refP, refQ, numberPSide, numberQSide, dbCoeffsP, dbCoeffsQ,
                             tc);
  }
}

/**
- Deblocking for the luminance component with strong or weak filter
.
\param piSrc           pointer to picture data
\param iOffset         offset value for picture data
\param tc              tc value
\param sw              decision strong/weak filter
\param bPartPNoFilter  indicator to disable filtering on partP
\param bPartQNoFilter  indicator to disable filtering on partQ
\param iThrCut         threshold value for weak filter decision
\param bFilterSecondP  decision weak filter/no filter for partP
\param bFilterSecondQ  decision weak filter/no filter for partQ
\param bitDepthLuma    luma bit depth
*/
template <typename TSrcDst>
void xPelFilterLumaCorePel(TSrcDst* piSrc, const ptrdiff_t iOffset, const int tc, const bool sw, const int iThrCut,
                           const bool bFilterSecondP, const bool bFilterSecondQ, const ClpRng& clpRng) {
  const int m0 = piSrc[-4 * iOffset];
  const int m1 = piSrc[-3 * iOffset];
  const int m2 = piSrc[-2 * iOffset];
  const int m3 = piSrc[-1 * iOffset];
  const int m4 = piSrc[0];
  const int m5 = piSrc[1 * iOffset];
  const int m6 = piSrc[2 * iOffset];
  const int m7 = piSrc[3 * iOffset];

  const char tc3[3] = {3, 2, 1};

  if (sw) {
    piSrc[-3 * iOffset] = Clip3(m1 - tc3[2] * tc, m1 + tc3[2] * tc, (2 * m0 + 3 * m1 + m2 + m3 + m4 + 4) >> 3);
    piSrc[-2 * iOffset] = Clip3(m2 - tc3[1] * tc, m2 + tc3[1] * tc, (m1 + m2 + m3 + m4 + 2) >> 2);
    piSrc[-1 * iOffset] = Clip3(m3 - tc3[0] * tc, m3 + tc3[0] * tc, (m1 + 2 * m2 + 2 * m3 + 2 * m4 + m5 + 4) >> 3);
    piSrc[0] = Clip3(m4 - tc3[0] * tc, m4 + tc3[0] * tc, (m2 + 2 * m3 + 2 * m4 + 2 * m5 + m6 + 4) >> 3);
    piSrc[1 * iOffset] = Clip3(m5 - tc3[1] * tc, m5 + tc3[1] * tc, (m3 + m4 + m5 + m6 + 2) >> 2);
    piSrc[2 * iOffset] = Clip3(m6 - tc3[2] * tc, m6 + tc3[2] * tc, (m3 + m4 + m5 + 3 * m6 + 2 * m7 + 4) >> 3);
  } else {
    /* Weak filter */
    int delta = (9 * (m4 - m3) - 3 * (m5 - m2) + 8) >> 4;

    if (abs(delta) < iThrCut) {
      delta = Clip3(-tc, tc, delta);
      const int tc2 = tc >> 1;

      piSrc[-iOffset * 1] = ClipPel(m3 + delta, clpRng);
      if (bFilterSecondP) {
        const int delta1 = Clip3(-tc2, tc2, ((((m1 + m3 + 1) >> 1) - m2 + delta) >> 1));
        piSrc[-iOffset * 2] = ClipPel(m2 + delta1, clpRng);
      }

      piSrc[0] = ClipPel(m4 - delta, clpRng);
      if (bFilterSecondQ) {
        const int delta2 = Clip3(-tc2, tc2, ((((m6 + m4 + 1) >> 1) - m5 - delta) >> 1));
        piSrc[iOffset] = ClipPel(m5 + delta2, clpRng);
      }
    }
  }
}

template <typename TSrcDst>
static inline void xPelFilterLumaCore(Pel* _piSrc, const ptrdiff_t step, const ptrdiff_t offset, const int tc,
                                      const bool sw, const int iThrCut, const bool bFilterSecondP,
                                      const bool bFilterSecondQ, const ClpRng& clpRng) {
  TSrcDst* piSrc = reinterpret_cast<TSrcDst*>(_piSrc);
  for (int i = 0; i < DEBLOCK_SMALLEST_BLOCK / 2; i++) {
    xPelFilterLumaCorePel<TSrcDst>(piSrc + step * i, offset, tc, sw, iThrCut, bFilterSecondP, bFilterSecondQ, clpRng);
  }
}

template <typename TSrcDst>
static inline void xEdgeFilterLumaImpCore(Pel* _piSrc, ptrdiff_t srcStep, ptrdiff_t offset, bool sidePisLarge,
                                          bool sideQisLarge, int iBeta, int iTc, int maxFilterLengthP,
                                          int maxFilterLengthQ, int iSideThreshold, int iThrCut, const ClpRng& clpRng) {
  TSrcDst* piSrc = reinterpret_cast<TSrcDst*>(_piSrc);
  const TSrcDst* piSrc0 = piSrc;
  const TSrcDst* piSrc3 = piSrc + 3 * srcStep;

  const int dp0 = xCalcDP<false, TSrcDst>(piSrc0, offset);
  const int dq0 = xCalcDQ<TSrcDst>(piSrc0, offset);
  const int dp3 = xCalcDP<false, TSrcDst>(piSrc3, offset);
  const int dq3 = xCalcDQ<TSrcDst>(piSrc3, offset);
  const int d0 = dp0 + dq0;
  const int d3 = dp3 + dq3;

  if (sidePisLarge || sideQisLarge) {
    const ptrdiff_t off3 = 3 * offset;
    const int dp0L = sidePisLarge ? ((dp0 + xCalcDP<false, TSrcDst>(piSrc0 - off3, offset) + 1) >> 1) : dp0;
    const int dq0L = sideQisLarge ? ((dq0 + xCalcDQ<TSrcDst>(piSrc0 + off3, offset) + 1) >> 1) : dq0;
    const int dp3L = sidePisLarge ? ((dp3 + xCalcDP<false, TSrcDst>(piSrc3 - off3, offset) + 1) >> 1) : dp3;
    const int dq3L = sideQisLarge ? ((dq3 + xCalcDQ<TSrcDst>(piSrc3 + off3, offset) + 1) >> 1) : dq3;

    const int d0L = dp0L + dq0L;
    const int d3L = dp3L + dq3L;

    const int dL = d0L + d3L;

    if (dL < iBeta) {
      // adjust decision so that it is not read beyond p5 is maxFilterLengthP is 5 and q5 if maxFilterLengthQ is 5
      const bool swL = xUseStrongFiltering<TSrcDst>(piSrc0, offset, 2 * d0L, iBeta, iTc, sidePisLarge, sideQisLarge,
                                                    maxFilterLengthP, maxFilterLengthQ) &&
                       xUseStrongFiltering<TSrcDst>(piSrc3, offset, 2 * d3L, iBeta, iTc, sidePisLarge, sideQisLarge,
                                                    maxFilterLengthP, maxFilterLengthQ);
      if (swL) {
        xFilteringPandQCore<TSrcDst>(reinterpret_cast<Pel*>(piSrc), srcStep, offset,
                                     sidePisLarge ? maxFilterLengthP : 3, sideQisLarge ? maxFilterLengthQ : 3, iTc);

        return;
      }
    }
  }

  // if( dL >= iBet || !swL )
  {
    const int dp = dp0 + dp3;
    const int dq = dq0 + dq3;
    const int d = d0 + d3;

    if (d < iBeta) {
      bool bFilterP = false;
      bool bFilterQ = false;

      if (maxFilterLengthP > 1 && maxFilterLengthQ > 1) {
        bFilterP = (dp < iSideThreshold);
        bFilterQ = (dq < iSideThreshold);
      }

      bool sw = false;

      if (maxFilterLengthP > 2 && maxFilterLengthQ > 2) {
        sw = xUseStrongFiltering<TSrcDst>(piSrc0, offset, 2 * d0, iBeta, iTc) &&
             xUseStrongFiltering<TSrcDst>(piSrc3, offset, 2 * d3, iBeta, iTc);
      }

      xPelFilterLumaCore<TSrcDst>(reinterpret_cast<Pel*>(piSrc), srcStep, offset, iTc, sw, iThrCut, bFilterP, bFilterQ,
                                  clpRng);
    }
  }
}

/**
- Deblocking of one line/column for the chrominance component
.
\param piSrc           pointer to picture data
\param iOffset         offset value for picture data
\param tc              tc value
\param bPartPNoFilter  indicator to disable filtering on partP
\param bPartQNoFilter  indicator to disable filtering on partQ
\param bitDepthChroma  chroma bit depth
*/
template <typename T>
static inline void xPelFilterChroma(T* piSrc, const ptrdiff_t iOffset, const int tc, const bool sw,
                                    const ClpRng& clpRng, const bool largeBoundary, const bool isChromaHorCTBBoundary) {
  int delta;

  const T m0 = piSrc[-iOffset * 4];
  const T m1 = piSrc[-iOffset * 3];
  const T m2 = piSrc[-iOffset * 2];
  const T m3 = piSrc[-iOffset * 1];
  const T m4 = piSrc[0];
  const T m5 = piSrc[iOffset * 1];
  const T m6 = piSrc[iOffset * 2];
  const T m7 = piSrc[iOffset * 3];

  if (sw) {
    if (isChromaHorCTBBoundary) {
      piSrc[-iOffset * 1] = Clip3(m3 - tc, m3 + tc, ((3 * m2 + 2 * m3 + m4 + m5 + m6 + 4) >> 3));      // p0
      piSrc[0] = Clip3(m4 - tc, m4 + tc, ((2 * m2 + m3 + 2 * m4 + m5 + m6 + m7 + 4) >> 3));            // q0
      piSrc[iOffset * 1] = Clip3(m5 - tc, m5 + tc, ((m2 + m3 + m4 + 2 * m5 + m6 + 2 * m7 + 4) >> 3));  // q1
      piSrc[iOffset * 2] = Clip3(m6 - tc, m6 + tc, ((m3 + m4 + m5 + 2 * m6 + 3 * m7 + 4) >> 3));       // q2
    } else {
      piSrc[-iOffset * 3] = Clip3(m1 - tc, m1 + tc, ((3 * m0 + 2 * m1 + m2 + m3 + m4 + 4) >> 3));        // p2
      piSrc[-iOffset * 2] = Clip3(m2 - tc, m2 + tc, ((2 * m0 + m1 + 2 * m2 + m3 + m4 + m5 + 4) >> 3));   // p1
      piSrc[-iOffset * 1] = Clip3(m3 - tc, m3 + tc, ((m0 + m1 + m2 + 2 * m3 + m4 + m5 + m6 + 4) >> 3));  // p0
      piSrc[0] = Clip3(m4 - tc, m4 + tc, ((m1 + m2 + m3 + 2 * m4 + m5 + m6 + m7 + 4) >> 3));             // q0
      piSrc[iOffset * 1] = Clip3(m5 - tc, m5 + tc, ((m2 + m3 + m4 + 2 * m5 + m6 + 2 * m7 + 4) >> 3));    // q1
      piSrc[iOffset * 2] = Clip3(m6 - tc, m6 + tc, ((m3 + m4 + m5 + 2 * m6 + 3 * m7 + 4) >> 3));         // q2
    }
  } else {
    delta = Clip3(-tc, tc, ((((m4 - m3) << 2) + m2 - m5 + 4) >> 3));

    piSrc[-iOffset] = ClipPel(m3 + delta, clpRng);
    piSrc[0] = ClipPel(m4 - delta, clpRng);
  }
}
// ====================================================================================================================
// Constructor / destructor / create / destroy
// ====================================================================================================================

void LoopFilter::initLoopFilter(int bitDepth) {
#define INIT_DB_FUNC(type)                           \
  xCalBsSameRef = xCalBsSameRefCore;                 \
  xEdgeFilterLumaImp = xEdgeFilterLumaImpCore<type>;

#if ADAPTIVE_BIT_DEPTH
  if (bitDepth <= 8) {
    INIT_DB_FUNC(Pel8bit);
  } else {
    INIT_DB_FUNC(Pel);
  }
#else
  INIT_DB_FUNC(Pel);
#endif

#undef INIT_DB_FUNC
}

LoopFilter::LoopFilter() {}

LoopFilter::~LoopFilter() {}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

/**
 - call deblocking function for every CU
 .
 \param  pcPic   picture class (Pic) pointer
 */

void LoopFilter::loopFilterCTU(CodingStructure& cs, const ChannelType chType, const int ctuCol, const int ctuLine,
                               const DeblockEdgeDir edgeDir, LoopFilterParam* loopFilterParamVerEdge) {
  PROFILER_SCOPE_AND_STAGE(1, g_timeProfiler, P_DBFILTER);
  const PreCalcValues& pcv = *cs.pcv;

  const int ly = ctuLine * pcv.maxCUHeight;

  if (ly >= pcv.lumaHeight) {
    return;
  }

  const UnitArea ctuArea = clipArea(
      UnitArea(pcv.chrFormat, Area(ctuCol << pcv.maxCUWidthLog2, ly, pcv.maxCUWidth, pcv.maxCUHeight)), *cs.picture);

  if (edgeDir == NUM_EDGE_DIR || edgeDir == EDGE_VER) {
    xDeblockCtuArea<EDGE_VER>(cs, ctuArea, chType, loopFilterParamVerEdge);
  }
  if (edgeDir == NUM_EDGE_DIR || edgeDir == EDGE_HOR) {
    xDeblockCtuArea<EDGE_HOR>(cs, ctuArea, chType, nullptr);
  }

#if ENABLE_TRACING
  for (unsigned x = 0; x < pcv.widthInCtus; x++) {
    const UnitArea ctuArea(pcv.chrFormat, Area(x << pcv.maxCUWidthLog2, ly, pcv.maxCUWidth, pcv.maxCUWidth));
    DTRACE(g_trace_ctx, D_CRC, "CTU %d %d", ctuArea.Y().x, ctuArea.Y().y);
    DTRACE_CRC(g_trace_ctx, D_CRC, cs, cs.picture->getRecoBuf(clipArea(ctuArea, *cs.picture)), &ctuArea.Y());
  }
#endif
}
// ====================================================================================================================
// Protected member functions
// ====================================================================================================================

/**
 Deblocking filter process in CU-based (the same function as conventional's)

 \param cu               the CU to be deblocked
 \param edgeDir          the direction of the edge in block boundary (horizontal/vertical), which is added newly
*/
template <DeblockEdgeDir edgeDir>
void LoopFilter::xDeblockCtuArea(CodingStructure& cs, const UnitArea& area, const ChannelType chType,
                                 LoopFilterParam* loopFilterParamVerEdge) {
  if (cs.getCtuCuPtrData(cs.ctuRsAddr(area.lumaPos(), CH_L)).cuPtr[0][0]->slice->getDeblockingFilterDisable()) {
    return;
  }

  const PreCalcValues& pcv = *cs.pcv;

  bool doLuma = chType == MAX_NUM_CHANNEL_TYPE || isLuma(chType);
#if JVET_Q0438_MONOCHROME_BUGFIXES
  bool doChroma = (chType == MAX_NUM_CHANNEL_TYPE || isChroma(chType)) && pcv.chrFormat != CHROMA_400 &&
                  area.blocks[COMPONENT_Cb].valid();
#else
  bool doChroma = (chType == MAX_NUM_CHANNEL_TYPE || isChroma(chType)) && pcv.chrFormat != CHROMA_400;
#endif
  static constexpr int incx = 4;
  static constexpr int incy = 4;

  const int csx = getChannelTypeScaleX(CH_C, pcv.chrFormat);
  const int csy = getChannelTypeScaleY(CH_C, pcv.chrFormat);

  LoopFilterParam const* lfpPtr;
  if (edgeDir == EDGE_HOR) {
    lfpPtr = cs.getLFPMapPtrHorEdge(cs.ctuRsAddr(area.Y().pos(), CH_L));
  } else {
    lfpPtr = loopFilterParamVerEdge;
  }

  ptrdiff_t lfpStride = cs.getLFPMapStride();

  const Position lumaPos = area.lumaPos();
  const Position chrmPos = doChroma ? area.chromaPos() : lumaPos;

#if ADAPTIVE_BIT_DEPTH
  if (cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) <= 8) {
    for (int dy = 0, cdy = 0; dy < area.lheight(); dy += incy, cdy += (incy >> csy)) {
      LoopFilterParam const* lineLfpPtr = lfpPtr;

      const int dyInCtu = (chrmPos.y + cdy) & (pcv.maxCUHeightMask >> csy);

      for (int dx = 0, cdx = 0; dx < area.lwidth(); dx += incx, cdx += (incx >> csx)) {
        const uint8_t bs = lineLfpPtr->bs;

        if (!doLuma || BsGet(bs, COMPONENT_Y) == 0) {
        } else {
          xEdgeFilterLuma<edgeDir, Pel8bit>(cs, {lumaPos.x + dx, lumaPos.y + dy}, *lineLfpPtr);
        }

        const int dxInCtu = (chrmPos.x + cdx) & (pcv.maxCUWidthMask >> csx);

        if (doChroma && (((edgeDir == EDGE_VER ? dxInCtu : dyInCtu) & (DEBLOCK_SMALLEST_BLOCK - 1)) == 0) &&
            (bs & 0x3c)) {
          xEdgeFilterChroma<edgeDir, Pel8bit>(cs, {chrmPos.x + cdx, chrmPos.y + cdy}, *lineLfpPtr);
        }

        OFFSETX(lineLfpPtr, lfpStride, 1);
      }

      OFFSETY(lfpPtr, lfpStride, 1);
    }
  } else {
    for (int dy = 0, cdy = 0; dy < area.lheight(); dy += incy, cdy += (incy >> csy)) {
      LoopFilterParam const* lineLfpPtr = lfpPtr;

      const int dyInCtu = (chrmPos.y + cdy) & (pcv.maxCUHeightMask >> csy);

      for (int dx = 0, cdx = 0; dx < area.lwidth(); dx += incx, cdx += (incx >> csx)) {
        const uint8_t bs = lineLfpPtr->bs;

        if (!doLuma || BsGet(bs, COMPONENT_Y) == 0) {
        } else {
          xEdgeFilterLuma<edgeDir, Pel>(cs, {lumaPos.x + dx, lumaPos.y + dy}, *lineLfpPtr);
        }

        const int dxInCtu = (chrmPos.x + cdx) & (pcv.maxCUWidthMask >> csx);

        if (doChroma && (((edgeDir == EDGE_VER ? dxInCtu : dyInCtu) & (DEBLOCK_SMALLEST_BLOCK - 1)) == 0) &&
            (bs & 0x3c)) {
          xEdgeFilterChroma<edgeDir, Pel>(cs, {chrmPos.x + cdx, chrmPos.y + cdy}, *lineLfpPtr);
        }

        OFFSETX(lineLfpPtr, lfpStride, 1);
      }

      OFFSETY(lfpPtr, lfpStride, 1);
    }
  }
#else
  for (int dy = 0, cdy = 0; dy < area.lheight(); dy += incy, cdy += (incy >> csy)) {
    LoopFilterParam const* lineLfpPtr = lfpPtr;

    const int dyInCtu = (chrmPos.y + cdy) & (pcv.maxCUHeightMask >> csy);

    for (int dx = 0, cdx = 0; dx < area.lwidth(); dx += incx, cdx += (incx >> csx)) {
      const uint8_t bs = lineLfpPtr->bs;

      if (!doLuma || BsGet(bs, COMPONENT_Y) == 0) {
      } else {
        xEdgeFilterLuma<edgeDir, Pel>(cs, {lumaPos.x + dx, lumaPos.y + dy}, *lineLfpPtr);
      }

      const int dxInCtu = (chrmPos.x + cdx) & (pcv.maxCUWidthMask >> csx);

      if (doChroma && (((edgeDir == EDGE_VER ? dxInCtu : dyInCtu) & (DEBLOCK_SMALLEST_BLOCK - 1)) == 0) &&
          (bs & 0x3c)) {
        xEdgeFilterChroma<edgeDir, Pel>(cs, {chrmPos.x + cdx, chrmPos.y + cdy}, *lineLfpPtr);
      }

      OFFSETX(lineLfpPtr, lfpStride, 1);
    }

    OFFSETY(lfpPtr, lfpStride, 1);
  }
#endif
}

void LoopFilter::calcFilterStrengthsHorEdge(const CodingUnit& cu) {
  if (cu.slice->getDeblockingFilterDisable()) {
    return;
  }

  const PreCalcValues& pcv = *cu.cs->pcv;
  const Area area = cu.blocks[cu.chType()];

  bool horEdgeFilter = false;
  int numHorVirBndry = 0;
  int horVirBndryPos[] = {0, 0, 0};

  const uint8_t channelScaleX = getChannelTypeScaleX(cu.chType(), cu.chromaFormat);
  const uint8_t channelScaleY = getChannelTypeScaleY(cu.chType(), cu.chromaFormat);

  const bool isCuCrossedByVirtualBoundariesHorEdge = isCrossedByVirtualBoundariesHorEdge(
      cu.cs->picHeader, area.y << channelScaleY, area.height << channelScaleY, numHorVirBndry, horVirBndryPos);

  if (isCuCrossedByVirtualBoundariesHorEdge) {
    CHECKD(numHorVirBndry >= (int)(sizeof(horVirBndryPos) / sizeof(horVirBndryPos[0])), "Too many virtual boundaries");
  }

  static constexpr int subBlockSize = 8;
  const PredictionUnit& currPU = cu;
  const Area& areaPu = area;
  LFCUParam stLFCUParam{xGetLoopfilterParam(cu)};
  const UnitScale scaling = cu.cs->getScaling(UnitScale::LF_PARAM_MAP, cu.chType());
  CtuData& ctuData = *cu.ctuData;
  // for SUBPU ATMVP and Affine, more PU deblocking needs to be found, for ISP the chroma block will be deferred to the
  // last luma block, so the processing order is different. For all other cases the boundary strenght can be directly
  // obtained in the TU loop.
  const bool refineBs =
      (currPU.mergeFlag() && currPU.mergeType() == MRG_TYPE_SUBPU_ATMVP) || currPU.affineFlag() || cu.ispMode();

  const int maskBlkY = ~((1 << scaling.posy) - 1);

  const bool pqCuSameCtuHor = (area.y & (pcv.maxCUHeightMask >> getChannelTypeScaleY(cu.chType(), pcv.chrFormat))) > 0;

  for (auto currTU = &cu.firstTU; currTU; currTU = currTU->next) {
    const Area& areaTu = currTU->blocks[currTU->chType];
    horEdgeFilter = (areaTu.y & maskBlkY) == area.y ? stLFCUParam.topEdge : true;

    const bool pqSameCtuHor = (areaTu.y & maskBlkY) == area.y ? pqCuSameCtuHor : true;

    if (isCuCrossedByVirtualBoundariesHorEdge) {
      xDeriveEdgefilterParam((areaTu.y & maskBlkY) << channelScaleY, numHorVirBndry, horVirBndryPos, horEdgeFilter);
    }

    xSetMaxFilterLengthPQFromTransformSizes<EDGE_HOR>(cu, *currTU, horEdgeFilter, !refineBs, ctuData, pqSameCtuHor,
                                                      nullptr);
  }

  if (!refineBs) return;

  if ((currPU.mergeFlag() && currPU.mergeType() == MRG_TYPE_SUBPU_ATMVP) || currPU.affineFlag()) {
    CHECK(cu.chType() != CH_L, "This path is only valid for single tree blocks!");

    for (int off = subBlockSize; off < areaPu.height; off += subBlockSize) {
      const Area mvBlockH(cu.Y().x, cu.Y().y + off, cu.Y().width, subBlockSize);
      horEdgeFilter = true;
      if (isCuCrossedByVirtualBoundariesHorEdge) {
        xDeriveEdgefilterParam(mvBlockH.y, numHorVirBndry, horVirBndryPos, horEdgeFilter);
      }

      xSetEdgeFilterInsidePu<EDGE_HOR>(cu, mvBlockH, horEdgeFilter, ctuData, nullptr);
    }

    xSetMaxFilterLengthPQForCodingSubBlocks<EDGE_HOR>(cu, ctuData, nullptr);
  }

  const unsigned uiPelsInPartX = pcv.minCUWidth >> channelScaleX;
  const unsigned uiPelsInPartY = pcv.minCUHeight >> channelScaleY;
  const ptrdiff_t lfpPos = cu.cs->inCtuPos(area.pos(), cu.chType());

  const CodingUnit* cuP = cu.above;
  const ChannelType chType = cu.chType();

  {
    LoopFilterParam* lfpPtrH = cu.cs->getLFPMapPtrHorEdge(cu.cs->ctuRsAddr(area.pos(), cu.chType()));
    ptrdiff_t lfpStride = cu.cs->getLFPMapStride();
    lfpPtrH += lfpPos;

    for (int y = 0; y < area.height; y += uiPelsInPartY) {
      LoopFilterParam* lineLfpPtrH = lfpPtrH;

      for (int x = 0; x < area.width; x += uiPelsInPartX) {
        cuP = (y || (cuP && cuP->blocks[chType].x + cuP->blocks[chType].width > area.x + x))
                  ? cuP
                  : cu.cs->getCU(Position{area.x + x, area.y - 1}, chType);

        if (lineLfpPtrH->filterEdge(cu.chType()))
          xGetBoundaryStrengthSingle<EDGE_HOR>(*lineLfpPtrH, cu, Position{area.x + x, area.y + y}, y ? cu : *cuP,
                                               ctuData, y ? true : pqCuSameCtuHor);

        lineLfpPtrH->bs &= ~BsSet(3, MAX_NUM_COMPONENT);

        INCX(lineLfpPtrH, lfpStride);
      }

      INCY(lfpPtrH, lfpStride);
    }
  }
}

void LoopFilter::calcFilterStrengthsVerEdge(const CodingUnit& cu, LoopFilterParam* loopFilterParamVerEdge) {
  if (cu.slice->getDeblockingFilterDisable()) {
    return;
  }

  const PreCalcValues& pcv = *cu.cs->pcv;
  const Area area = cu.blocks[cu.chType()];

  bool verEdgeFilter = false;
  int numVerVirBndry = 0;
  int verVirBndryPos[] = {0, 0, 0};

  const uint8_t channelScaleX = getChannelTypeScaleX(cu.chType(), cu.chromaFormat);
  const uint8_t channelScaleY = getChannelTypeScaleY(cu.chType(), cu.chromaFormat);

  const bool isCuCrossedByVirtualBoundariesVerEdge = isCrossedByVirtualBoundariesVerEdge(
      cu.cs->picHeader, area.x << channelScaleX, area.width << channelScaleX, numVerVirBndry, verVirBndryPos);

  if (isCuCrossedByVirtualBoundariesVerEdge) {
    CHECKD(numVerVirBndry >= (int)(sizeof(verVirBndryPos) / sizeof(verVirBndryPos[0])), "Too many virtual boundaries");
  }

  static constexpr int subBlockSize = 8;
  const PredictionUnit& currPU = cu;
  const Area& areaPu = area;
  LFCUParam stLFCUParam{xGetLoopfilterParam(cu)};
  const UnitScale scaling = cu.cs->getScaling(UnitScale::LF_PARAM_MAP, cu.chType());
  CtuData& ctuData = *cu.ctuData;
  // for SUBPU ATMVP and Affine, more PU deblocking needs to be found, for ISP the chroma block will be deferred to the
  // last luma block, so the processing order is different. For all other cases the boundary strenght can be directly
  // obtained in the TU loop.
  const bool refineBs =
      (currPU.mergeFlag() && currPU.mergeType() == MRG_TYPE_SUBPU_ATMVP) || currPU.affineFlag() || cu.ispMode();

  const int maskBlkX = ~((1 << scaling.posx) - 1);

  const bool pqCuSameCtuVer = (area.x & (pcv.maxCUWidthMask >> getChannelTypeScaleX(cu.chType(), pcv.chrFormat))) > 0;

  for (auto currTU = &cu.firstTU; currTU; currTU = currTU->next) {
    const Area& areaTu = currTU->blocks[currTU->chType];

    verEdgeFilter = (areaTu.x & maskBlkX) == area.x ? stLFCUParam.leftEdge : true;

    const bool pqSameCtuVer = (areaTu.x & maskBlkX) == area.x ? pqCuSameCtuVer : true;

    if (isCuCrossedByVirtualBoundariesVerEdge) {
      xDeriveEdgefilterParam((areaTu.x & maskBlkX) << channelScaleX, numVerVirBndry, verVirBndryPos, verEdgeFilter);
    }

    xSetMaxFilterLengthPQFromTransformSizes<EDGE_VER>(cu, *currTU, verEdgeFilter, !refineBs, ctuData, pqSameCtuVer,
                                                      loopFilterParamVerEdge);
  }

  if (!refineBs) return;

  if ((currPU.mergeFlag() && currPU.mergeType() == MRG_TYPE_SUBPU_ATMVP) || currPU.affineFlag()) {
    CHECK(cu.chType() != CH_L, "This path is only valid for single tree blocks!");

    for (int off = subBlockSize; off < areaPu.width; off += subBlockSize) {
      const Area mvBlockV(cu.Y().x + off, cu.Y().y, subBlockSize, cu.Y().height);
      verEdgeFilter = true;
      if (isCuCrossedByVirtualBoundariesVerEdge) {
        xDeriveEdgefilterParam(mvBlockV.x, numVerVirBndry, verVirBndryPos, verEdgeFilter);
      }

      xSetEdgeFilterInsidePu<EDGE_VER>(cu, mvBlockV, verEdgeFilter, ctuData, loopFilterParamVerEdge);
    }

    xSetMaxFilterLengthPQForCodingSubBlocks<EDGE_VER>(cu, ctuData, loopFilterParamVerEdge);
  }

  const unsigned uiPelsInPartX = pcv.minCUWidth >> channelScaleX;
  const unsigned uiPelsInPartY = pcv.minCUHeight >> channelScaleY;
  const ptrdiff_t lfpPos = cu.cs->inCtuPos(area.pos(), cu.chType());

  const CodingUnit* cuP = cu.left;
  const ChannelType chType = cu.chType();

  {
    LoopFilterParam* lfpPtrV = loopFilterParamVerEdge;
    ptrdiff_t lfpStride = cu.cs->getLFPMapStride();
    lfpPtrV += lfpPos;

    for (int y = 0; y < area.height; y += uiPelsInPartY) {
      LoopFilterParam* lineLfpPtrV = lfpPtrV;

      cuP = cuP && cuP->blocks[chType].y + cuP->blocks[chType].height > area.y + y
                ? cuP
                : cu.cs->getCU(Position{area.x - 1, area.y + y}, chType);

      for (int x = 0; x < area.width; x += uiPelsInPartX) {
        if (lineLfpPtrV->filterEdge(cu.chType()))
          xGetBoundaryStrengthSingle<EDGE_VER>(*lineLfpPtrV, cu, Position{area.x + x, area.y + y}, x ? cu : *cuP,
                                               ctuData, x ? true : pqCuSameCtuVer);

        lineLfpPtrV->bs &= ~BsSet(3, MAX_NUM_COMPONENT);

        INCX(lineLfpPtrV, lfpStride);
      }

      INCY(lfpPtrV, lfpStride);
    }
  }
}

void LoopFilter::calcFilterStrengthsVerEdgeCTU(CodingStructure& cs, const UnitArea& ctuArea,
                                               LoopFilterParam* loopFilterParamVerEdge) {
  for (auto& currCU : cs.traverseCUs(ctuArea)) {
    calcFilterStrengthsVerEdge(currCU, loopFilterParamVerEdge);
  }
}

inline bool LoopFilter::isCrossedByVirtualBoundariesHorEdge(const PicHeader* picHeader, int y, int height,
                                                            int& numHorVirBndry, int horVirBndryPos[]) {
  numHorVirBndry = 0;
  if (!picHeader->getVirtualBoundariesPresentFlag()) {
    return false;
  }

  for (int i = 0; i < picHeader->getNumHorVirtualBoundaries(); i++) {
    if (y <= picHeader->getVirtualBoundariesPosY(i) && picHeader->getVirtualBoundariesPosY(i) < y + height) {
      horVirBndryPos[numHorVirBndry++] = picHeader->getVirtualBoundariesPosY(i);
    }
  }

  return numHorVirBndry > 0;
}

inline bool LoopFilter::isCrossedByVirtualBoundariesVerEdge(const PicHeader* picHeader, int x, int width,
                                                            int& numVerVirBndry, int verVirBndryPos[]) {
  numVerVirBndry = 0;
  if (!picHeader->getVirtualBoundariesPresentFlag()) {
    return false;
  }

  for (int i = 0; i < picHeader->getNumVerVirtualBoundaries(); i++) {
    if (x <= picHeader->getVirtualBoundariesPosX(i) && picHeader->getVirtualBoundariesPosX(i) < x + width) {
      verVirBndryPos[numVerVirBndry++] = picHeader->getVirtualBoundariesPosX(i);
    }
  }

  return numVerVirBndry > 0;
}

inline void LoopFilter::xDeriveEdgefilterParam(const int pos, const int numVirBndry, const int virBndryPos[],
                                               bool& edgeFilter) {
  for (int i = 0; i < numVirBndry; i++) {
    if (virBndryPos[i] == pos) {
      edgeFilter = false;
      break;
    }
  }
}

template <DeblockEdgeDir edgeDir>
inline PosType parlPos(const Position&);
template <DeblockEdgeDir edgeDir>
inline PosType perpPos(const Position&);

template <DeblockEdgeDir edgeDir>
inline SizeType parlSize(const Size&);
template <DeblockEdgeDir edgeDir>
inline SizeType perpSize(const Size&);

template <>
inline PosType parlPos<EDGE_HOR>(const Position& pos) {
  return pos.x;
}
template <>
inline PosType perpPos<EDGE_HOR>(const Position& pos) {
  return pos.y;
}
template <>
inline PosType parlPos<EDGE_VER>(const Position& pos) {
  return pos.y;
}
template <>
inline PosType perpPos<EDGE_VER>(const Position& pos) {
  return pos.x;
}

template <>
inline SizeType parlSize<EDGE_HOR>(const Size& size) {
  return size.width;
}
template <>
inline SizeType perpSize<EDGE_HOR>(const Size& size) {
  return size.height;
}
template <>
inline SizeType parlSize<EDGE_VER>(const Size& size) {
  return size.height;
}
template <>
inline SizeType perpSize<EDGE_VER>(const Size& size) {
  return size.width;
}

template <DeblockEdgeDir edgeDir>
void LoopFilter::xSetMaxFilterLengthPQForCodingSubBlocks(const CodingUnit& cu, CtuData& ctuData,
                                                         LoopFilterParam* loopFilterParamVerEdge) {
  static constexpr int subBlockSize = 8;
  const int minCUWidth = cu.cs->pcv->minCUWidth;
  const int minCUHeight = cu.cs->pcv->minCUHeight;
  const int xInc = edgeDir ? minCUWidth : subBlockSize;
  const int yInc = edgeDir ? subBlockSize : minCUHeight;

  LoopFilterParam* lfpPtrL;

  if (edgeDir == EDGE_HOR) {
    lfpPtrL = cu.cs->getLFPMapPtrHorEdge(cu.lumaPos(), CH_L);
  } else {
    lfpPtrL = loopFilterParamVerEdge + cu.cs->inCtuPos(cu.lumaPos(), CH_L);
  }
  ptrdiff_t lfpStride = cu.cs->getLFPMapStride();
  const UnitScale scaling = cu.cs->getScaling(UnitScale::LF_PARAM_MAP, CHANNEL_TYPE_LUMA);

  for (int y = 0; y < cu.Y().height; y += yInc) {
    LoopFilterParam* lfpPtr = lfpPtrL;

    for (int x = 0; x < cu.Y().width; x += xInc) {
      const PosType perpVal = edgeDir ? y : x;

      if (lfpPtr->edge) {
        lfpPtr->sideMaxFiltLengthQ = std::min<int>(lfpPtr->sideMaxFiltLengthQ, 5);

        if (perpVal > 0) {
          lfpPtr->sideMaxFiltLengthP = std::min<int>(lfpPtr->sideMaxFiltLengthP, 5);
        }
      } else if (perpVal > 0 && (GET_OFFSET(lfpPtr, lfpStride, -(1 - edgeDir), -edgeDir)->edge ||
                                 (perpVal + 4 >= perpSize<edgeDir>(cu.Y())) ||
                                 GET_OFFSET(lfpPtr, lfpStride, (1 - edgeDir), edgeDir)->edge)) {
        lfpPtr->sideMaxFiltLengthQ = 1;
        lfpPtr->sideMaxFiltLengthP = 1;
      } else if (perpVal > 0 && ((edgeDir ? (y == 8) : (x == 8)) ||
                                 GET_OFFSET(lfpPtr, lfpStride, -2 * (1 - edgeDir), -2 * edgeDir)->edge ||
                                 (perpVal + 8 >= perpSize<edgeDir>(cu.Y())) ||
                                 GET_OFFSET(lfpPtr, lfpStride, 2 * (1 - edgeDir), 2 * edgeDir)->edge)) {
        lfpPtr->sideMaxFiltLengthQ = 2;
        lfpPtr->sideMaxFiltLengthP = 2;
      } else {
        lfpPtr->sideMaxFiltLengthQ = 3;
        lfpPtr->sideMaxFiltLengthP = 3;
      }

      OFFSETX(lfpPtr, lfpStride, scaling.scaleHor(xInc));
    }

    OFFSETY(lfpPtrL, lfpStride, scaling.scaleVer(yInc));
  }
}

static inline TransformUnit const* getTU(const CodingUnit& cu, const Position& pos, const ChannelType chType) {
  const TransformUnit* ptu = &cu.firstTU;

  while (!(ptu->blocks[chType].x + ptu->blocks[chType].width > pos.x &&
           ptu->blocks[chType].y + ptu->blocks[chType].height > pos.y)) {
    ptu = ptu->next;
  }

  return ptu;
}

template <DeblockEdgeDir edgeDir>
void LoopFilter::xSetMaxFilterLengthPQFromTransformSizes(const CodingUnit& cu, const TransformUnit& currTU,
                                                         const bool bValue, bool deriveBdStrngt, CtuData& ctuData,
                                                         bool pqSameCtu, LoopFilterParam* loopFilterParamVerEdge) {
  const PreCalcValues& pcv = *cu.cs->pcv;

  ChannelType start = CH_L;
  ChannelType end = CH_C;

  const bool dt = CU::isSepTree(cu);

  if (dt) {
    if (cu.chType() == CH_L) {
      end = CH_L;
    } else {
      start = CH_C;
    }
  }

  if (start != end && (!isChromaEnabled(pcv.chrFormat) || !currTU.Cb().valid())) {
    end = CH_L;
  }

  const int csx = (start != end) ? getChannelTypeScaleX(CH_C, pcv.chrFormat) : 0;
  const int csy = (start != end) ? getChannelTypeScaleY(CH_C, pcv.chrFormat) : 0;
  const Area& area = currTU.blocks[end];

  for (int ct = start; ct <= end; ct++) {
    const ChannelType ch = (ChannelType)ct;
    const bool vld = isLuma(ch) ? currTU.Y().valid() : currTU.Cb().valid();
    if (vld && perpPos<edgeDir>(currTU.blocks[ch]) != 0) {
      LoopFilterParam* lfpPtr;
      if (edgeDir == EDGE_HOR) {
        lfpPtr = cu.cs->getLFPMapPtrHorEdge(currTU.blocks[ct], ch);
      } else {
        lfpPtr = loopFilterParamVerEdge + cu.cs->inCtuPos(currTU.blocks[ct], ch);
      }

      ptrdiff_t lfpStride = cu.cs->getLFPMapStride();

      const int inc = edgeDir ? pcv.minCUWidth >> getChannelTypeScaleX(ch, cu.chromaFormat)
                              : pcv.minCUHeight >> getChannelTypeScaleY(ch, cu.chromaFormat);

      const CodingUnit* cuNeigh = edgeDir ? cu.above : cu.left;
      const CodingUnit* cuP =
          (cuNeigh && perpPos<edgeDir>(currTU.blocks[ch]) == perpPos<edgeDir>(cu.blocks[ch])) ? cuNeigh : &cu;
      const CodingUnit* cuPfstCh =
          (cuNeigh && perpPos<edgeDir>(currTU.blocks[start]) == perpPos<edgeDir>(cu.blocks[start])) ? cuNeigh : &cu;
      const int incFst = edgeDir ? pcv.minCUWidth >> getChannelTypeScaleX(ChannelType(start), cu.chromaFormat)
                                 : pcv.minCUHeight >> getChannelTypeScaleY(ChannelType(start), cu.chromaFormat);

      if (cuP == &cu && (edgeDir ? cu.above : cu.left) == nullptr &&
          (edgeDir ? cu.blocks[ch].pos().y : cu.blocks[ch].pos().x) >
              0)  // TODO: check for !pps.getLoopFilterAcrossSlicesEnabledFlag() ||
                  // !pps.getLoopFilterAcrossTilesEnabledFlag()
      {
        const Position posP{currTU.blocks[ch].x - (1 - edgeDir), currTU.blocks[ch].y - edgeDir};
        const Position posPfst{currTU.blocks[start].x - (1 - edgeDir), currTU.blocks[start].y - edgeDir};
        cuP = cuP->blocks[ch].contains(posP) ? cuP : cu.cs->getCU(posP, ch);
        cuPfstCh = cuPfstCh->blocks[start].contains(posPfst) ? cuPfstCh : cu.cs->getCU(posPfst, start);
      }
      bool bSameCUTUSize = perpPos<edgeDir>(currTU.blocks[ch]) == perpPos<edgeDir>(cu.blocks[ch]);
      if (ct == end && deriveBdStrngt) {
        if (start != end) {
          for (int d = 0, dFst = 0; d < parlSize<edgeDir>(currTU.blocks[ch]);) {
            const Position posQ{currTU.blocks[ch].x + edgeDir * d, currTU.blocks[ch].y + (1 - edgeDir) * d};
            const Position posP = posQ.offset(-(1 - edgeDir), -edgeDir);
            const int sizeQSide = perpSize<edgeDir>(currTU.blocks[ch]);
            cuP = parlPos<edgeDir>(cuP->blocks[ch]) + parlSize<edgeDir>(cuP->blocks[ch]) > parlPos<edgeDir>(posP)
                      ? cuP
                      : cu.cs->getCU(posP, ch);
            const Position posPfst = currTU.blocks[start].offset(edgeDir ? dFst : -1, edgeDir ? -1 : dFst);
            cuPfstCh = parlPos<edgeDir>(cuPfstCh->blocks[start]) + parlSize<edgeDir>(cuPfstCh->blocks[start]) >
                               parlPos<edgeDir>(posPfst)
                           ? cuPfstCh
                           : cu.cs->getCU(posPfst, ChannelType(start));
            const TransformUnit& tuP = cuP->firstTU.next == nullptr ? cuP->firstTU : *getTU(*cuP, posP, ch);
            const int sizePSide = perpSize<edgeDir>(tuP.blocks[ch]);
            LoopFilterParam& lfp = *lfpPtr;
            lfp.setFilterCMFL((sizeQSide >= 8 && sizePSide >= 8) ? 1 : 0);
            if (bValue)
              xGetBoundaryStrengthSingle<edgeDir>(
                  lfp, cu, Position((area.x + edgeDir * d) << csx, (area.y + (1 - edgeDir) * d) << csy), *cuPfstCh,
                  ctuData, pqSameCtu);
            lfp.bs &= ~BsSet(3, MAX_NUM_COMPONENT);
            if (!CU::isIntra(cu) && !CU::isIntra(*cuP) && cuP == cuPfstCh && cu.geoFlag() == false &&
                cuP->geoFlag() == false) {
              int distance =
                  parlPos<edgeDir>(tuP.blocks[ch]) + parlSize<edgeDir>(tuP.blocks[ch]) - parlPos<edgeDir>(posQ);
              if (distance > parlSize<edgeDir>(currTU.blocks[ch]) - d)
                distance = parlSize<edgeDir>(currTU.blocks[ch]) - d;  // range check
              if (distance > inc && !cuP->affineFlag() &&
                  !(cuP->mergeFlag() && cuP->mergeType() == MRG_TYPE_SUBPU_ATMVP)) {
                LoopFilterParam param = lfp;
                for (int j = inc; j < distance; j += inc, d += inc, dFst += incFst) {
                  OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
                  *lfpPtr = param;
                }
              }
            }
            OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
            d += inc;
            dFst += incFst;
          }
        } else {
          // old logical
          for (int d = 0, dFst = 0; d < parlSize<edgeDir>(currTU.blocks[ch]); d += inc, dFst += incFst) {
            const Position posQ{currTU.blocks[ch].x + edgeDir * d, currTU.blocks[ch].y + (1 - edgeDir) * d};
            const Position posP = posQ.offset(-(1 - edgeDir), -edgeDir);
            const Position posPfst = currTU.blocks[start].offset(edgeDir ? dFst : -1, edgeDir ? -1 : dFst);
            const int sizeQSide = perpSize<edgeDir>(currTU.blocks[ch]);
            cuP = parlPos<edgeDir>(cuP->blocks[ch]) + parlSize<edgeDir>(cuP->blocks[ch]) > parlPos<edgeDir>(posP)
                      ? cuP
                      : cu.cs->getCU(posP, ch);
            cuPfstCh = start != end
                           ? parlPos<edgeDir>(cuPfstCh->blocks[start]) + parlSize<edgeDir>(cuPfstCh->blocks[start]) >
                                     parlPos<edgeDir>(posPfst)
                                 ? cuPfstCh
                                 : cu.cs->getCU(posPfst, ChannelType(start))
                           : cuP;
            const TransformUnit& tuP = cuP->firstTU.next == nullptr ? cuP->firstTU : *getTU(*cuP, posP, ch);
            const int sizePSide = perpSize<edgeDir>(tuP.blocks[ch]);
            LoopFilterParam& lfp = *lfpPtr;
            lfp.setFilterEdge(cu.chType(), bValue);
            if ((lfp.bs || bSameCUTUSize) && bValue) {
              lfp.bs |= BsSet(3, MAX_NUM_COMPONENT);
            } else {
              lfp.bs |= BsSet(1, MAX_NUM_COMPONENT);
            }
            if (ch == CHANNEL_TYPE_LUMA) {
              bool smallBlock = (sizePSide <= 4) || (sizeQSide <= 4);
              if (smallBlock) {
                lfp.sideMaxFiltLengthP = 1;
                lfp.sideMaxFiltLengthQ = 1;
              } else {
                lfp.sideMaxFiltLengthP = (sizePSide >= 32) ? (cuP->affineFlag() ? 5 : 7) : 3;
                lfp.sideMaxFiltLengthQ = (sizeQSide >= 32) ? 7 : 3;
              }
              lfp.edge = 1;
            } else {
              lfp.setFilterCMFL((sizeQSide >= 8 && sizePSide >= 8) ? 1 : 0);
            }

            if (bValue)
              xGetBoundaryStrengthSingle<edgeDir>(
                  lfp, cu, Position((area.x + edgeDir * d) << csx, (area.y + (1 - edgeDir) * d) << csy), *cuPfstCh,
                  ctuData, pqSameCtu);
            lfp.bs &= ~BsSet(3, MAX_NUM_COMPONENT);
            OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
          }
        }

      } else {
        for (int d = 0; d < parlSize<edgeDir>(currTU.blocks[ch]);) {
          const Position posQ{currTU.blocks[ch].x + edgeDir * d, currTU.blocks[ch].y + (1 - edgeDir) * d};
          const Position posP = posQ.offset(-(1 - edgeDir), -edgeDir);
          const int sizeQSide = perpSize<edgeDir>(currTU.blocks[ch]);
          cuP = parlPos<edgeDir>(cuP->blocks[ch]) + parlSize<edgeDir>(cuP->blocks[ch]) > parlPos<edgeDir>(posP)
                    ? cuP
                    : cu.cs->getCU(posP, ch);
          const TransformUnit& tuP = cuP->firstTU.next == nullptr ? cuP->firstTU : *getTU(*cuP, posP, ch);
          const int sizePSide = perpSize<edgeDir>(tuP.blocks[ch]);
          int distance = parlPos<edgeDir>(tuP.blocks[ch]) + parlSize<edgeDir>(tuP.blocks[ch]) - parlPos<edgeDir>(posQ);
          if (distance > parlSize<edgeDir>(currTU.blocks[ch]) - d)
            distance = parlSize<edgeDir>(currTU.blocks[ch]) - d;  // range check
          if (ch == CHANNEL_TYPE_LUMA) {
            LoopFilterParam& lfp = *lfpPtr;
            lfp.setFilterEdge(cu.chType(), bValue);
            if ((lfp.bs || bSameCUTUSize) && bValue) {
              lfp.bs |= BsSet(3, MAX_NUM_COMPONENT);
            } else {
              lfp.bs |= BsSet(1, MAX_NUM_COMPONENT);
            }
            bool smallBlock = (sizePSide <= 4) || (sizeQSide <= 4);
            if (smallBlock) {
              lfp.sideMaxFiltLengthP = 1;
              lfp.sideMaxFiltLengthQ = 1;
            } else {
              lfp.sideMaxFiltLengthP = (sizePSide >= 32) ? (cuP->affineFlag() ? 5 : 7) : 3;
              lfp.sideMaxFiltLengthQ = (sizeQSide >= 32) ? 7 : 3;
            }
            lfp.edge = 1;
            if (distance > inc) {
              LoopFilterParam tmp = lfp;
              for (int j = inc; j < distance; j += inc, d += inc) {
                OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
                *lfpPtr = tmp;
              }
            }
            d += inc;
            OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
          } else {
            for (int j = 0; j < distance; j += inc, d += inc) {
              LoopFilterParam& lfp = *lfpPtr;
              if (ct == start) {
                lfp.setFilterEdge(cu.chType(), bValue);
                if ((lfp.bs || bSameCUTUSize) && bValue) {
                  lfp.bs |= BsSet(3, MAX_NUM_COMPONENT);
                } else {
                  lfp.bs |= BsSet(1, MAX_NUM_COMPONENT);
                }
              }
              lfp.setFilterCMFL((sizeQSide >= 8 && sizePSide >= 8) ? 1 : 0);
              OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
            }
          }
        }
      }
    }
  }
}

template <DeblockEdgeDir edgeDir>
void LoopFilter::xSetEdgeFilterInsidePu(const CodingUnit& cu, const Area& area, const bool bValue, CtuData& ctuData,
                                        LoopFilterParam* loopFilterParamVerEdge) {
  const PreCalcValues& pcv = *cu.cs->pcv;
  LoopFilterParam* lfpPtr;
  if (edgeDir == EDGE_HOR) {
    lfpPtr = cu.cs->getLFPMapPtrHorEdge(area.pos(), cu.chType());
  } else {
    lfpPtr = loopFilterParamVerEdge + cu.cs->inCtuPos(area.pos(), cu.chType());
  }
  ptrdiff_t lfpStride = cu.cs->getLFPMapStride();

  const int inc = edgeDir ? pcv.minCUWidth >> getChannelTypeScaleX(cu.chType(), cu.chromaFormat)
                          : pcv.minCUHeight >> getChannelTypeScaleY(cu.chType(), cu.chromaFormat);
  if (bValue) {
    for (int d = 0; d < parlSize<edgeDir>(area); d += inc) {
      lfpPtr->setFilterEdge(cu.chType(), 1);
      if (lfpPtr->bs) {
        lfpPtr->bs |= BsSet(3, MAX_NUM_COMPONENT);
      }
      OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
    }
  } else {
    for (int d = 0; d < parlSize<edgeDir>(area); d += inc) {
      lfpPtr->setFilterEdge(cu.chType(), 0);
      OFFSET(lfpPtr, lfpStride, edgeDir, (1 - edgeDir));
    }
  }
}

LFCUParam LoopFilter::xGetLoopfilterParam(const CodingUnit& cu) {
  const Position pos = cu.blocks[cu.chType()].pos();

  const PPS& pps = *cu.cs->pps;
  const SPS& sps = *cu.cs->sps;

  const CodingUnit& cuLeft =
      (pos.x > 0 && cu.left == nullptr) ? *cu.cs->getCU(pos.offset(-1, 0), cu.chType()) : *cu.left;
  const CodingUnit& cuAbove =
      (pos.y > 0 && cu.above == nullptr) ? *cu.cs->getCU(pos.offset(0, -1), cu.chType()) : *cu.above;

  const bool loopFilterAcrossSubPicEnabledFlagLeft =
      !sps.getSubPicInfoPresentFlag() || (pps.getSubPicFromCU(cu).getloopFilterAcrossSubPicEnabledFlag() &&
                                          pps.getSubPicFromCU(cuLeft).getloopFilterAcrossSubPicEnabledFlag());
  const bool loopFilterAcrossSubPicEnabledFlagTop =
      !sps.getSubPicInfoPresentFlag() || (pps.getSubPicFromCU(cu).getloopFilterAcrossSubPicEnabledFlag() &&
                                          pps.getSubPicFromCU(cuAbove).getloopFilterAcrossSubPicEnabledFlag());

  LFCUParam stLFCUParam;  ///< status structure
  stLFCUParam.leftEdge =
      (0 < pos.x) && isAvailable(cu, cuLeft, !pps.getLoopFilterAcrossSlicesEnabledFlag(),
                                 !pps.getLoopFilterAcrossTilesEnabledFlag(), !loopFilterAcrossSubPicEnabledFlagLeft);
  stLFCUParam.topEdge =
      (0 < pos.y) && isAvailable(cu, cuAbove, !pps.getLoopFilterAcrossSlicesEnabledFlag(),
                                 !pps.getLoopFilterAcrossTilesEnabledFlag(), !loopFilterAcrossSubPicEnabledFlagTop);
  return stLFCUParam;
}

template <DeblockEdgeDir edgeDir>
void LoopFilter::xGetBoundaryStrengthSingle(LoopFilterParam& lfp, const CodingUnit& cuQ, const Position& localPos,
                                            const CodingUnit& cuP, CtuData& ctuData, bool pqSameCtu) {
  const Slice& sliceQ = *cuQ.slice;
  const ChannelType chType = cuQ.chType();
  const Position& cuPos = cuQ.blocks[chType].pos();
  const Position& posQ = localPos;
  const Position posP{posQ.x - !edgeDir, posQ.y - edgeDir};

  const TransformUnit& tuQ = cuQ.firstTU.next == nullptr ? cuQ.firstTU : *getTU(cuQ, posQ, chType);
  const TransformUnit& tuP =
      cuP.firstTU.next == nullptr
          ? cuP.firstTU
          : *getTU(cuP, posP, chType);  // TODO: check this: based on chType of the current cu, because cuQ.chType and
                                        // cuP.chType are not the same when local dual-tree is applied

  const bool hasLuma = cuQ.Y().valid();
  const bool hasChroma = isChromaEnabled(cuQ.chromaFormat) && cuQ.Cb().valid();

  bool cuPcIsIntra = false;
  int chrmBS = 2;

  if (hasLuma) {
    lfp.qp[0] = (cuQ.qp + cuP.qp + 1) >> 1;
  }

  if (hasChroma) {
    const int qpBdOffset2 = cuQ.cs->sps->getQpBDOffset(CH_C) << 1;
    const bool isPQDiffCh = !chType && cuP.treeType() != TREE_D;
    const TransformUnit& tuQc = cuQ.ispMode() ? *cuQ.lastTU : tuQ;
    const Position posPc = isPQDiffCh ? recalcPosition(cuQ.chromaFormat, chType, CH_C, posP) : Position();
    const CodingUnit& cuPc = isPQDiffCh ? *cuQ.cs->getCU(posPc, CH_C) : cuP;
    const TransformUnit& tuPc = isPQDiffCh ? *getTU(cuPc, posPc, CH_C) : (cuP.ispMode() ? *cuP.lastTU : tuP);

    cuPcIsIntra = CU::isIntra(cuPc);

    lfp.qp[1] = ((tuPc.chromaQp[0] + tuQc.chromaQp[0] - qpBdOffset2 + 1) >> 1);
    lfp.qp[2] = ((tuPc.chromaQp[1] + tuQc.chromaQp[1] - qpBdOffset2 + 1) >> 1);

    if (cuPcIsIntra) {
      chrmBS = (MODE_INTRA == cuPc.predMode() && cuPc.bdpcmModeChroma()) &&
                       (MODE_INTRA == cuQ.predMode() && cuQ.bdpcmModeChroma())
                   ? 0
                   : 2;
    }
  }

  const int bsMask = (hasLuma ? BsSet(3, COMPONENT_Y) : 0) | BsSet(3, MAX_NUM_COMPONENT) |
                     (hasChroma ? BsSet(3, COMPONENT_Cb) : 0) | (hasChroma ? BsSet(3, COMPONENT_Cr) : 0);

  //-- Set BS for Intra MB : BS = 4 or 3
  if (MODE_INTRA == cuP.predMode() || MODE_INTRA == cuQ.predMode()) {
    const int edgeIdx = (perpPos<edgeDir>(localPos) - perpPos<edgeDir>(cuPos)) / 4;

    int bsY = cuP.bdpcmMode() && cuQ.bdpcmMode() ? 0 : 2;

    if (cuQ.ispMode() && edgeIdx) {
      lfp.bs |= BsSet(bsY, COMPONENT_Y) & bsMask;
    } else {
      lfp.bs |= (BsSet(bsY, COMPONENT_Y) + BsSet(chrmBS, COMPONENT_Cb) + BsSet(chrmBS, COMPONENT_Cr)) & bsMask;
    }

    return;
  } else if (cuPcIsIntra) {
    lfp.bs |= (BsSet(chrmBS, COMPONENT_Cb) + BsSet(chrmBS, COMPONENT_Cr));
  }

  if ((lfp.bs & bsMask) && (cuP.ciipFlag() || cuQ.ciipFlag())) {
    lfp.bs |= (BsSet(2, COMPONENT_Y) + BsSet(2, COMPONENT_Cb) + BsSet(2, COMPONENT_Cr)) & bsMask;

    return;
  }

  unsigned tmpBs = 0;
  //-- Set BS for not Intra MB : BS = 2 or 1 or 0
  // Y
  if (lfp.bs & bsMask) {
    int cbfSum = tuQ.cbf | tuP.cbf;
    tmpBs += BsSet(cbfSum & 1, COMPONENT_Y);
    if (!(MODE_INTRA != cuP.predMode() && MODE_INTRA != cuQ.predMode() && cuPcIsIntra)) {
      bool jointChr = tuQ.jointCbCr || tuP.jointCbCr;
      cbfSum >>= 1;
      tmpBs += BsSet((cbfSum & 1) | jointChr, COMPONENT_Cb);
      cbfSum >>= 1;
      tmpBs += BsSet((cbfSum & 1) | jointChr, COMPONENT_Cr);
    }
  }

  if (BsGet(tmpBs, COMPONENT_Y) == 1) {
    lfp.bs |= tmpBs & bsMask;

    return;
  }

  if (cuP.ciipFlag() || cuQ.ciipFlag()) {
    lfp.bs |= 1 & bsMask;

    return;
  }

  if (!hasLuma) {
    lfp.bs |= tmpBs & bsMask;
    return;
  }

  if (BsGet(lfp.bs, MAX_NUM_COMPONENT) != 0 && BsGet(lfp.bs, MAX_NUM_COMPONENT) != 3) {
    lfp.bs |= tmpBs & bsMask;
    return;
  }

  if (hasChroma) {
    lfp.bs |= tmpBs & bsMask;
  }

  if (cuP.predMode() != cuQ.predMode() && hasLuma) {
    lfp.bs |= 1 & bsMask;
    return;
  }

  const ptrdiff_t pqDiff = edgeDir ? int(cuQ.cs->getLFPMapStride()) : 1;
  const ptrdiff_t inCtuPosQ = cuQ.cs->inCtuPos(posQ, chType);

  // and now the pred
  const MotionInfo& miQ = ctuData.motion[inCtuPosQ];
  const MotionInfo& miP = !pqSameCtu ? cuP.getMotionInfo(posP) : ctuData.motion[inCtuPosQ - pqDiff];
  const Slice& sliceP = *cuP.slice;

  static constexpr int nThreshold = (1 << MV_FRACTIONAL_BITS_INTERNAL) >> 1;

  if (sliceQ.isInterB() || sliceP.isInterB()) {
    //    const Picture* piRefP0 = CU::isIBC(cuP)       ? sliceP.getPic()
    //                             : 0 <= miP.refIdx[0] ? sliceP.getRefPic(REF_PIC_LIST_0, miP.refIdx[0])
    //                                                  : nullptr;
    //    const Picture* piRefP1 = CU::isIBC(cuP)       ? nullptr
    //                             : 0 <= miP.refIdx[1] ? sliceP.getRefPic(REF_PIC_LIST_1, miP.refIdx[1])
    //                                                  : nullptr;
    //    const Picture* piRefQ0 = CU::isIBC(cuQ)       ? sliceQ.getPic()
    //                             : 0 <= miQ.refIdx[0] ? sliceQ.getRefPic(REF_PIC_LIST_0, miQ.refIdx[0])
    //                                                  : nullptr;
    //    const Picture* piRefQ1 = CU::isIBC(cuQ)       ? nullptr
    //                             : 0 <= miQ.refIdx[1] ? sliceQ.getRefPic(REF_PIC_LIST_1, miQ.refIdx[1])
    //                                                  : nullptr;
    const Picture* piRefP0 = 0 <= miP.refIdx[0] ? sliceP.getRefPic(REF_PIC_LIST_0, miP.refIdx[0]) : nullptr;
    const Picture* piRefP1 = 0 <= miP.refIdx[1] ? sliceP.getRefPic(REF_PIC_LIST_1, miP.refIdx[1]) : nullptr;
    const Picture* piRefQ0 = 0 <= miQ.refIdx[0] ? sliceQ.getRefPic(REF_PIC_LIST_0, miQ.refIdx[0]) : nullptr;
    const Picture* piRefQ1 = 0 <= miQ.refIdx[1] ? sliceQ.getRefPic(REF_PIC_LIST_1, miQ.refIdx[1]) : nullptr;

    unsigned uiBs = 0;

    // th can be optimized
    if ((piRefP0 == piRefQ0 && piRefP1 == piRefQ1) || (piRefP0 == piRefQ1 && piRefP1 == piRefQ0)) {
      uiBs = xCalBsSameRef(miP, miQ, piRefP0, piRefP1, piRefQ0, piRefQ1);
    } else  // for all different Ref_Idx
    {
      uiBs = 1;
    }

    lfp.bs |= (uiBs + tmpBs) & bsMask;

    return;
  }

  // pcSlice->isInterP()
  // CHECK( CU::isInter( cuP ) && 0 > miP.refIdx[0], "Invalid reference picture list index" );
  // CHECK( CU::isInter( cuP ) && 0 > miQ.refIdx[0], "Invalid reference picture list index" );

  const Picture* piRefP0 = (CU::isIBC(cuP) ? sliceP.getPic() : sliceP.getRefPic(REF_PIC_LIST_0, miP.refIdx[0]));
  const Picture* piRefQ0 = (CU::isIBC(cuQ) ? sliceQ.getPic() : sliceQ.getRefPic(REF_PIC_LIST_0, miQ.refIdx[0]));

  if (piRefP0 != piRefQ0) {
    lfp.bs |= (tmpBs + 1) & bsMask;

    return;
  }

  Mv mvP0 = miP.mv[0];
  Mv mvQ0 = miQ.mv[0];

  lfp.bs |= (((abs(mvQ0.getHor() - mvP0.getHor()) >= nThreshold) || (abs(mvQ0.getVer() - mvP0.getVer()) >= nThreshold))
                 ? (tmpBs + 1)
                 : tmpBs) &
            bsMask;
}

#if LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET
template <typename T>
void LoopFilter::deriveLADFShift(const T* src, const ptrdiff_t stride, int& shift, const DeblockEdgeDir edgeDir,
                                 const SPS sps) {
  int32_t lumaLevel = 0;
  shift = sps.getLadfQpOffset(0);

  if (edgeDir == EDGE_VER) {
    lumaLevel = (src[0] + src[3 * stride] + src[-1] + src[3 * stride - 1]) >> 2;
  } else  // (edgeDir == EDGE_HOR)
  {
    lumaLevel = (src[0] + src[3] + src[-stride] + src[-stride + 3]) >> 2;
  }

  for (int k = 1; k < sps.getLadfNumIntervals(); k++) {
    const int th = sps.getLadfIntervalLowerBound(k);
    if (lumaLevel > th) {
      shift = sps.getLadfQpOffset(k);
    } else {
      break;
    }
  }
}
#endif

template <DeblockEdgeDir edgeDir, typename T>
void LoopFilter::xEdgeFilterLuma(CodingStructure& cs, const Position& pos, const LoopFilterParam& lfp) {
#if ADAPTIVE_BIT_DEPTH
  int bytePerPixel = cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) <= 8 ? 1 : 2;
#endif

  PelBuf picYuvRec = cs.getRecoBuf(COMPONENT_Y);
#if ADAPTIVE_BIT_DEPTH
  T* piSrc = reinterpret_cast<T*>(picYuvRec.bufAt(pos, bytePerPixel));
#else
  T* piSrc = picYuvRec.bufAt(pos);
#endif
  const ptrdiff_t iStride = picYuvRec.stride;
  const SPS& sps = *cs.sps;
  const Slice& slice = *cs.m_ctuCuPtr[cs.ctuRsAddr(pos, CH_L)].cuPtr[0][0]->slice;
  const int bitDepthLuma = sps.getBitDepth(CHANNEL_TYPE_LUMA);
  const ClpRng& clpRng = slice.clpRng(COMPONENT_Y);

  const int betaOffsetDiv2 = slice.getDeblockingFilterBetaOffsetDiv2();
  const int tcOffsetDiv2 = slice.getDeblockingFilterTcOffsetDiv2();

  ptrdiff_t offset, srcStep;

  if (edgeDir == EDGE_VER) {
    offset = 1;
    srcStep = iStride;
  } else  // (edgeDir == EDGE_HOR)
  {
    offset = iStride;
    srcStep = 1;
  }

#if ENABLE_SIMD_OPT
  if (edgeDir == EDGE_VER) {
    _mm_prefetch((char*)&piSrc[0 * iStride - 4], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[1 * iStride - 4], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[2 * iStride - 4], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[3 * iStride - 4], _MM_HINT_T0);
  } else {
    _mm_prefetch((char*)&piSrc[(0 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(1 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(2 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(3 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(4 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(5 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(6 - 4) * iStride], _MM_HINT_T0);
    _mm_prefetch((char*)&piSrc[(7 - 4) * iStride], _MM_HINT_T0);
  }
#endif  // ENABLE_SIMD_OPT

  const unsigned uiBs = BsGet(lfp.bs, COMPONENT_Y);

  if (!uiBs) {
    return;
  }

  int iQP = lfp.qp[0];
#if LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET
  if (sps.getLadfEnabled()) {
    int iShift = 0;
    deriveLADFShift<T>(piSrc, iStride, iShift, edgeDir, sps);
    iQP += iShift;
  }

#endif
  const int maxFilterLengthP = lfp.sideMaxFiltLengthP;
  const int maxFilterLengthQ = lfp.sideMaxFiltLengthQ;

  bool sidePisLarge = maxFilterLengthP > 3;
  bool sideQisLarge = maxFilterLengthQ > 3;

  if (edgeDir == EDGE_HOR && (pos.y & cs.pcv->maxCUHeightMask) == 0) {
    sidePisLarge = false;
  }

  const int iIndexTC =
      Clip3(0, MAX_QP + DEFAULT_INTRA_TC_OFFSET, int(iQP + DEFAULT_INTRA_TC_OFFSET * (uiBs - 1) + (tcOffsetDiv2 << 1)));
  const int iIndexB = Clip3(0, MAX_QP, iQP + (betaOffsetDiv2 << 1));

  const int iTc = bitDepthLuma < 10 ? ((sm_tcTable[iIndexTC] + (1 << (9 - bitDepthLuma))) >> (10 - bitDepthLuma))
                                    : ((sm_tcTable[iIndexTC]) << (bitDepthLuma - 10));
  const int iBeta = sm_betaTable[iIndexB] << (bitDepthLuma - 8);
  const int iSideThreshold = (iBeta + (iBeta >> 1)) >> 3;
  const int iThrCut = iTc * 10;

  static constexpr bool bPartPNoFilter = false;
  static constexpr bool bPartQNoFilter = false;

  if (!bPartPNoFilter || !bPartQNoFilter) {
    xEdgeFilterLumaImp(reinterpret_cast<Pel*>(piSrc), srcStep, offset, sidePisLarge, sideQisLarge, iBeta, iTc,
                       maxFilterLengthP, maxFilterLengthQ, iSideThreshold, iThrCut, clpRng);
  }
}

template <DeblockEdgeDir edgeDir, typename T>
void LoopFilter::xEdgeFilterChroma(CodingStructure& cs, const Position& pos, const LoopFilterParam& lfp) {
  const PreCalcValues& pcv = *cs.pcv;

#if ADAPTIVE_BIT_DEPTH
  int bytePerPixel = cs.sps->getBitDepth(CHANNEL_TYPE_LUMA) <= 8 ? 1 : 2;
#endif

  const ChromaFormat nChromaFormat = pcv.chrFormat;
  const int csy = getChannelTypeScaleY(CH_C, nChromaFormat);
  const unsigned uiPelsInPartChromaH = pcv.minCUWidth >> getChannelTypeScaleX(CH_C, nChromaFormat);
  const unsigned uiPelsInPartChromaV = pcv.minCUHeight >> csy;

  PelBuf picYuvRecCb = cs.getRecoBuf(COMPONENT_Cb);
  PelBuf picYuvRecCr = cs.getRecoBuf(COMPONENT_Cr);
#if ADAPTIVE_BIT_DEPTH
  T* piSrcCb = reinterpret_cast<T*>(picYuvRecCb.bufAt(pos, bytePerPixel));
  T* piSrcCr = reinterpret_cast<T*>(picYuvRecCr.bufAt(pos, bytePerPixel));
#else
  T* piSrcCb = picYuvRecCb.bufAt(pos);
  T* piSrcCr = picYuvRecCr.bufAt(pos);
#endif
  const ptrdiff_t iStride = picYuvRecCb.stride;
  const SPS& sps = *cs.sps;
  const Slice& slice = *cs.m_ctuCuPtr[cs.ctuRsAddr(pos, CH_C)].cuPtr[1][0]->slice;
  const int bitDepthChroma = sps.getBitDepth(CHANNEL_TYPE_CHROMA);

  const int tcOffsetDiv2[2] = {slice.getDeblockingFilterCbTcOffsetDiv2(), slice.getDeblockingFilterCrTcOffsetDiv2()};
  const int betaOffsetDiv2[2] = {slice.getDeblockingFilterCbBetaOffsetDiv2(),
                                 slice.getDeblockingFilterCrBetaOffsetDiv2()};

  ptrdiff_t offset, srcStep;
  unsigned uiLoopLength;

  if (edgeDir == EDGE_VER) {
    offset = 1;
    srcStep = iStride;
    uiLoopLength = uiPelsInPartChromaV;
  } else {
    offset = iStride;
    srcStep = 1;
    uiLoopLength = uiPelsInPartChromaH;
  }

  unsigned bS[2];

  unsigned tmpBs = lfp.bs;
  bS[0] = BsGet(tmpBs, COMPONENT_Cb);
  bS[1] = BsGet(tmpBs, COMPONENT_Cr);

  if (bS[0] <= 0 && bS[1] <= 0) {
    return;
  }

  bool largeBoundary = lfp.filterCMFL();
  bool isChromaHorCTBBoundary = false;

  if (edgeDir == EDGE_HOR && (pos.y & (pcv.maxCUHeightMask >> csy)) == 0) {
    isChromaHorCTBBoundary = true;
  }

  static constexpr bool bPartPNoFilter = false;
  static constexpr bool bPartQNoFilter = false;

  if (!bPartPNoFilter || !bPartQNoFilter)
    for (unsigned chromaIdx = 0; chromaIdx < 2; chromaIdx++) {
      if (bS[chromaIdx] == 2 || (largeBoundary && bS[chromaIdx] == 1)) {
        const ClpRng& clpRng(cs.picture->slices[0]->clpRng(ComponentID(chromaIdx + 1)));

        int iQP = lfp.qp[chromaIdx + 1];

        const int iIndexTC =
            Clip3<int>(0, MAX_QP + DEFAULT_INTRA_TC_OFFSET,
                       iQP + DEFAULT_INTRA_TC_OFFSET * (bS[chromaIdx] - 1) + (tcOffsetDiv2[chromaIdx] << 1));
        const int iTc = bitDepthChroma < 10
                            ? ((sm_tcTable[iIndexTC] + (1 << (9 - bitDepthChroma))) >> (10 - bitDepthChroma))
                            : ((sm_tcTable[iIndexTC]) << (bitDepthChroma - 10));
        T* piSrcChroma = chromaIdx == 0 ? piSrcCb : piSrcCr;

        if (largeBoundary) {
          const int iBitdepthScale = 1 << (sps.getBitDepth(CHANNEL_TYPE_CHROMA) - 8);

          const int indexB = Clip3<int>(0, MAX_QP, iQP + (betaOffsetDiv2[chromaIdx] << 1));
          const int beta = sm_betaTable[indexB] * iBitdepthScale;

          const int dp0 =
              isChromaHorCTBBoundary ? xCalcDP<true, T>(piSrcChroma, offset) : xCalcDP<false, T>(piSrcChroma, offset);
          const int dq0 = xCalcDQ<T>(piSrcChroma, offset);
          const int subSamplingShift = (edgeDir == EDGE_VER) ? getChannelTypeScaleY(CH_C, nChromaFormat)
                                                             : getChannelTypeScaleX(CH_C, nChromaFormat);
          const int dp3 = isChromaHorCTBBoundary
                              ? ((subSamplingShift == 1) ? xCalcDP<true, T>(piSrcChroma + srcStep, offset)
                                                         : xCalcDP<true, T>(piSrcChroma + srcStep * 3, offset))
                              : ((subSamplingShift == 1) ? xCalcDP<false, T>(piSrcChroma + srcStep, offset)
                                                         : xCalcDP<false, T>(piSrcChroma + srcStep * 3, offset));
          const int dq3 = (subSamplingShift == 1) ? xCalcDQ<T>(piSrcChroma + srcStep, offset)
                                                  : xCalcDQ<T>(piSrcChroma + srcStep * 3, offset);

          const int d0 = dp0 + dq0;
          const int d3 = dp3 + dq3;
          const int d = d0 + d3;

          if (d < beta) {
            const bool sw =
                xUseStrongFiltering<T>(piSrcChroma, offset, 2 * d0, beta, iTc, false, false, 7, 7,
                                       isChromaHorCTBBoundary) &&
                xUseStrongFiltering<T>(piSrcChroma + ((subSamplingShift == 1) ? srcStep : srcStep * 3), offset, 2 * d3,
                                       beta, iTc, false, false, 7, 7, isChromaHorCTBBoundary);

            for (unsigned uiStep = 0; uiStep < uiLoopLength; uiStep++) {
              xPelFilterChroma<T>(piSrcChroma + srcStep * uiStep, offset, iTc, sw, clpRng, largeBoundary,
                                  isChromaHorCTBBoundary);
            }

            continue;
          }
        }

        for (unsigned uiStep = 0; uiStep < uiLoopLength; uiStep++) {
          xPelFilterChroma<T>(piSrcChroma + srcStep * uiStep, offset, iTc, false, clpRng, largeBoundary,
                              isChromaHorCTBBoundary);
        }
      }
    }
}

//! \}
