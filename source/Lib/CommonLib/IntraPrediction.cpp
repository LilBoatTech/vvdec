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

/** \file     Prediction.cpp
    \brief    prediction class
*/

#define DONT_UNDEF_SIZE_AWARE_PER_EL_OP

#include "IntraPrediction.h"

#include "Unit.h"
#include "UnitTools.h"

#include "Buffer.h"

#include "dtrace_next.h"
#include "Rom.h"

#include <memory.h>
#include <array>

#include "CommonLib/InterpolationFilter.h"
#include "CommonLib/TimeProfiler.h"

//! \ingroup CommonLib
//! \{

// ====================================================================================================================
// Tables
// ====================================================================================================================

const uint8_t IntraPrediction::m_aucIntraFilter[MAX_NUM_CHANNEL_TYPE][MAX_INTRA_FILTER_DEPTHS] = {{
                                                                                                      // Luma
                                                                                                      24,  //   1xn
                                                                                                      24,  //   2xn
                                                                                                      24,  //   4xn
                                                                                                      14,  //   8xn
                                                                                                      2,   //  16xn
                                                                                                      0,   //  32xn
                                                                                                      0,   //  64xn
                                                                                                      0,   // 128xn
                                                                                                  },
                                                                                                  {
                                                                                                      // Chroma
                                                                                                      40,  //   1xn
                                                                                                      40,  //   2xn
                                                                                                      40,  //   4xn
                                                                                                      28,  //   8xn
                                                                                                      4,   //  16xn
                                                                                                      0,   //  32xn
                                                                                                      0,   //  64xn
                                                                                                      0,   // 128xn
                                                                                                  }};

const TFilterCoeff g_intraGaussFilter[32][4] = {
    {16, 32, 16, 0}, {16, 32, 16, 0}, {15, 31, 17, 1}, {15, 31, 17, 1}, {14, 30, 18, 2}, {14, 30, 18, 2},
    {13, 29, 19, 3}, {13, 29, 19, 3}, {12, 28, 20, 4}, {12, 28, 20, 4}, {11, 27, 21, 5}, {11, 27, 21, 5},
    {10, 26, 22, 6}, {10, 26, 22, 6}, {9, 25, 23, 7},  {9, 25, 23, 7},  {8, 24, 24, 8},  {8, 24, 24, 8},
    {7, 23, 25, 9},  {7, 23, 25, 9},  {6, 22, 26, 10}, {6, 22, 26, 10}, {5, 21, 27, 11}, {5, 21, 27, 11},
    {4, 20, 28, 12}, {4, 20, 28, 12}, {3, 19, 29, 13}, {3, 19, 29, 13}, {2, 18, 30, 14}, {2, 18, 30, 14},
    {1, 17, 31, 15}, {1, 17, 31, 15},
};

template <typename TSrcDst>
static void IntraHorVerPDPCCore(Pel *_pDstBuf, const ptrdiff_t dstStride, Pel *_refSide, const CPelBuf &pSrc,
                                const int width, const int height, int scale, const Pel *_refMain,
                                const ClpRng &clpRng) {
  TSrcDst *pDstBuf = reinterpret_cast<TSrcDst *>(_pDstBuf);
  TSrcDst *refSide = reinterpret_cast<TSrcDst *>(_refSide);
  const TSrcDst *refMain = reinterpret_cast<const TSrcDst *>(_refMain);
  const int lev[4] = {std::min(3, width), std::min(6, width), std::min(12, width), std::min(24, width)};

  const TSrcDst topLeft = (reinterpret_cast<const TSrcDst *>(pSrc.buf))[0];
  for (int y = 0; y < height; y++) {
    const int left = refSide[y + 1];
    TSrcDst *line = &pDstBuf[y * dstStride];
    for (int x = 0; x < lev[scale]; x++) {
      int wL = 32 >> std::min(31, ((x << 1) >> scale));
      *line++ = ClipPel((wL * (left - topLeft) + (refMain[x + 1] << 6) + 32) >> 6, clpRng);
    }
    memcpy(line, refMain + lev[scale] + 1, (width - lev[scale]) * sizeof(TSrcDst));
  }
}

template <typename TSrcDst>
static void IntraAnglePDPCCore(Pel *_pDsty, const ptrdiff_t dstStride, Pel *_refSide, const int width, const int height,
                               int scale, int invAngle) {
  TSrcDst *pDsty = reinterpret_cast<TSrcDst *>(_pDsty);
  TSrcDst *refSide = reinterpret_cast<TSrcDst *>(_refSide);
  for (int y = 0; y < height; y++, pDsty += dstStride) {
    int invAngleSum = 256;

    for (int x = 0; x < std::min(3 << scale, width); x++) {
      invAngleSum += invAngle;

      int wL = 32 >> (2 * x >> scale);
      int left = refSide[y + (invAngleSum >> 9) + 1];
      pDsty[x] = pDsty[x] + ((wL * (left - pDsty[x]) + 32) >> 6);
    }
  }
}

template <typename TSrcDst>
void GetLumaRecPixel420Core(const int width, const int height, const Pel *_pRecSrc0, const ptrdiff_t iRecStride,
                            Pel *_pDst0, const ptrdiff_t iDstStride) {
  const TSrcDst *pRecSrc0 = reinterpret_cast<const TSrcDst *>(_pRecSrc0);
  TSrcDst *pDst0 = reinterpret_cast<TSrcDst *>(_pDst0);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      pDst0[x + 0] =
          (pRecSrc0[((x + 0) << 1)] * 2 + pRecSrc0[((x + 0) << 1) + 1] * 1 + pRecSrc0[((x + 0) << 1) - 1] * 1 +
           pRecSrc0[((x + 0) << 1) + iRecStride] * 2 + pRecSrc0[((x + 0) << 1) + 1 + iRecStride] * 1 +
           pRecSrc0[((x + 0) << 1) - 1 + iRecStride] * 1 + 4) >>
          3;
    }
    pDst0 += iDstStride;
    pRecSrc0 += (iRecStride << 1);
  }
}

/** Function for deriving planar intra prediction. This function derives the prediction samples for planar mode (intra
 * coding).
 */

// NOTE: Bit-Limit - 24-bit source
template <typename TSrcDst>
void xPredIntraPlanarCore(const CPelBuf &pSrc, PelBuf &pDst, const SPS &sps) {
  // with some optimizations gcc gives spurious "-Wmaybe-uninitialized" warnings here
  GCC_WARNING_DISABLE_maybe_uninitialized

      const uint32_t width = pDst.width;
  const uint32_t height = pDst.height;
  const uint32_t log2W = getLog2(width);
  const uint32_t log2H = getLog2(height);
  int leftColumn[MAX_CU_SIZE + 1], topRow[MAX_CU_SIZE + 1], bottomRow[MAX_CU_SIZE], rightColumn[MAX_CU_SIZE];
  const uint32_t offset = 1 << (log2W + log2H);
  const ptrdiff_t srcStride = pSrc.stride;
  const TSrcDst *src = (const TSrcDst *)pSrc.buf;
  // Get left and above reference column and row
  for (int k = 0; k < width + 1; k++) {
    topRow[k] = src[(k + 1) + 0 * srcStride];
  }

  for (int k = 0; k < height + 1; k++) {
    leftColumn[k] = src[0 + (k + 1) * srcStride];
  }

  // Prepare intermediate variables used in interpolation
  int bottomLeft = leftColumn[height];
  int topRight = topRow[width];

  for (int k = 0; k < width; k++) {
    bottomRow[k] = bottomLeft - topRow[k];
    topRow[k] = topRow[k] << log2H;
  }

  for (int k = 0; k < height; k++) {
    rightColumn[k] = topRight - leftColumn[k];
    leftColumn[k] = leftColumn[k] << log2W;
  }

  const uint32_t finalShift = 1 + log2W + log2H;
  const ptrdiff_t stride = pDst.stride;
  TSrcDst *pred = reinterpret_cast<TSrcDst *>(pDst.buf);
  for (int y = 0; y < height; y++, pred += stride) {
    int horPred = leftColumn[y];

    for (int x = 0; x < width; x++) {
      horPred += rightColumn[y];
      topRow[x] += bottomRow[x];

      int vertPred = topRow[x];
      pred[x] = ((horPred << log2H) + (vertPred << log2W) + offset) >> finalShift;
    }
  }
  GCC_WARNING_RESET
}

template <typename TSrcDst>
void IntraPredSampleFilterCore(Pel *_ptrSrc, const ptrdiff_t srcStride, PelBuf &piPred, const uint32_t uiDirMode,
                               const ClpRng &clpRng) {
  if (uiDirMode == PLANAR_IDX || uiDirMode == DC_IDX) {
    const int iWidth = piPred.width;
    const int iHeight = piPred.height;
    TSrcDst *ptrSrc = reinterpret_cast<TSrcDst *>(_ptrSrc);
    TSrcDst *dstBuf = reinterpret_cast<TSrcDst *>(piPred.buf);
    int dstStride = static_cast<int>(piPred.stride);

    const int scale = ((getLog2(iWidth) - 2 + getLog2(iHeight) - 2 + 2) >> 2);
    CHECK(scale < 0 || scale > 31, "PDPC: scale < 0 || scale > 31");
    for (int y = 0; y < iHeight; y++) {
      const int wT = 32 >> std::min(31, ((y << 1) >> scale));
      //      const Pel left = srcBuf.at(y + 1, 1);
      const int left = ptrSrc[0 + (y + 1) * srcStride];
      for (int x = 0; x < iWidth; x++) {
        const int wL = 32 >> std::min(31, ((x << 1) >> scale));
        const int top = ptrSrc[(x + 1) + 0 * srcStride];
        const int val = dstBuf[y * dstStride + x];
        dstBuf[y * dstStride + x] = val + ((wL * (left - val) + wT * (top - val) + 32) >> 6);
      }
    }
  }
}

template <typename TSrcDst>
void IntraPredAngleLumaCore(Pel *_pDstBuf, const ptrdiff_t dstStride, Pel *_refMain, int width, int height,
                            int deltaPos, int intraPredAngle, const TFilterCoeff *ff, const bool useCubicFilter,
                            const ClpRng &clpRng) {
  TSrcDst *pDstBuf = reinterpret_cast<TSrcDst *>(_pDstBuf);
  TSrcDst *refMain = reinterpret_cast<TSrcDst *>(_refMain);
  for (int y = 0; y < height; y++) {
    const int deltaInt = deltaPos >> 5;
    const int deltaFract = deltaPos & (32 - 1);

    int p[4];

    int refMainIndex = deltaInt + 1;

    const TFilterCoeff *f = &ff[deltaFract << 2];

    for (int x = 0; x < width; x++, refMainIndex++) {
      p[0] = refMain[refMainIndex - 1];
      p[1] = refMain[refMainIndex];
      p[2] = refMain[refMainIndex + 1];
      p[3] = refMain[refMainIndex + 2];

      if (useCubicFilter) {  // only cubic filter has negative coefficients and requires clipping
        int temp = ((static_cast<int>(f[0] * p[0]) + static_cast<int>(f[1] * p[1]) + static_cast<int>(f[2] * p[2]) +
                     static_cast<int>(f[3] * p[3]) + 32) >>
                    6);
        pDstBuf[y * dstStride + x] = ClipPel(temp, clpRng);
      } else {
        pDstBuf[y * dstStride + x] = ((static_cast<int>(f[0] * p[0]) + static_cast<int>(f[1] * p[1]) +
                                       static_cast<int>(f[2] * p[2]) + static_cast<int>(f[3] * p[3]) + 32) >>
                                      6);
      }
    }
    deltaPos += intraPredAngle;
  }
}

template <typename TSrcDst>
void IntraPredAngleChromaCore(Pel *_pDstBuf, const ptrdiff_t dstStride, Pel *_pBorder, int width, int height,
                              int deltaPos, int intraPredAngle) {
  TSrcDst *pDstBuf = reinterpret_cast<TSrcDst *>(_pDstBuf);
  TSrcDst *pBorder = reinterpret_cast<TSrcDst *>(_pBorder);
  for (int y = 0; y < height; y++) {
    const int deltaInt = deltaPos >> 5;
    const int deltaFract = deltaPos & (32 - 1);

    // Do linear filtering
    const TSrcDst *pRM = pBorder + deltaInt + 1;
    int lastRefMainPel = *pRM++;

    for (int x = 0; x < width; pRM++, x++) {
      int thisRefMainPel = *pRM;
      pDstBuf[x + 0] = (((32 - deltaFract) * lastRefMainPel + deltaFract * thisRefMainPel + 16) >> 5);
      lastRefMainPel = thisRefMainPel;
    }
    deltaPos += intraPredAngle;
    pDstBuf += dstStride;
  }
}

// ====================================================================================================================
// Constructor / destructor / initialize
// ====================================================================================================================

IntraPrediction::IntraPrediction() : m_currChromaFormat(NUM_CHROMA_FORMAT) {
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++) {
    for (uint32_t buf = 0; buf < NUM_PRED_BUF; buf++) {
      m_piYuvExt[ch][buf] = nullptr;
    }
  }
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++) {
    for (uint32_t buf = 0; buf < 4; buf++) {
      m_yuvExt2[ch][buf] = nullptr;
    }
  }

  m_piTemp = nullptr;
  m_pMdlmTemp = nullptr;

#if ADAPTIVE_BIT_DEPTH
  m_bytePerPixel = 2;
#endif
}

IntraPrediction::~IntraPrediction() { destroy(); }

void IntraPrediction::initIntraPrediction(int bitDepth) {
#define INIT_INTRA_FUNC(type)                                \
  IntraPredSampleFilter8 = IntraPredSampleFilterCore<type>;  \
  IntraPredSampleFilter16 = IntraPredSampleFilterCore<type>; \
  IntraPredAngleCore = IntraPredAngleLumaCore<type>;         \
  IntraPredAngleCore4 = IntraPredAngleLumaCore<type>;        \
  IntraPredAngleCore8 = IntraPredAngleLumaCore<type>;        \
  IntraPredAngleChroma = IntraPredAngleChromaCore<type>;     \
  IntraPredAngleChroma4 = IntraPredAngleChromaCore<type>;    \
  IntraPredAngleChroma8 = IntraPredAngleChromaCore<type>;    \
  xPredIntraPlanar = xPredIntraPlanarCore<type>;             \
  IntraHorVerPDPC = IntraHorVerPDPCCore<type>;               \
  IntraAnglePDPC = IntraAnglePDPCCore<type>;                 \
  GetLumaRecPixel420 = GetLumaRecPixel420Core<type>;

#if ADAPTIVE_BIT_DEPTH
  m_bytePerPixel = bitDepth <= 8 ? 1 : 2;
  if (bitDepth <= 8) {
    INIT_INTRA_FUNC(Pel8bit);
  } else {
    INIT_INTRA_FUNC(Pel);
  }
#else
  INIT_INTRA_FUNC(Pel);
#endif

#undef INIT_INTRA_FUNC
}

void IntraPrediction::destroy() {
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++) {
    for (uint32_t buf = 0; buf < NUM_PRED_BUF; buf++) {
      delete[] m_piYuvExt[ch][buf];
      m_piYuvExt[ch][buf] = nullptr;
    }
  }
  for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++) {
    for (uint32_t buf = 0; buf < 4; buf++) {
      delete[] m_yuvExt2[ch][buf];
      m_yuvExt2[ch][buf] = nullptr;
    }
  }

  delete[] m_piTemp;
  m_piTemp = nullptr;
  delete[] m_pMdlmTemp;
  m_pMdlmTemp = nullptr;
}

void IntraPrediction::init(ChromaFormat chromaFormatIDC, const unsigned bitDepthY) {
  // if it has been initialised before, but the chroma format has changed, release the memory and start again.
  if (m_piYuvExt[COMPONENT_Y][PRED_BUF_UNFILTERED] != nullptr && m_currChromaFormat != chromaFormatIDC) {
    destroy();
  }

  if (m_yuvExt2[COMPONENT_Y][0] != nullptr && m_currChromaFormat != chromaFormatIDC) {
    destroy();
  }

  m_currChromaFormat = chromaFormatIDC;

  if (m_piYuvExt[COMPONENT_Y][PRED_BUF_UNFILTERED] ==
      nullptr)  // check if first is null (in which case, nothing initialised yet)
  {
    m_iYuvExtSize = (MAX_CU_SIZE * 2 + 1 + MAX_REF_LINE_IDX) * (MAX_CU_SIZE * 2 + 1 + MAX_REF_LINE_IDX);

    for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++) {
      for (uint32_t buf = 0; buf < NUM_PRED_BUF; buf++) {
        m_piYuvExt[ch][buf] = new Pel[m_iYuvExtSize];
      }
    }
  }

  if (m_yuvExt2[COMPONENT_Y][0] == nullptr)  // check if first is null (in which case, nothing initialised yet)
  {
    m_yuvExtSize2 = (MAX_CU_SIZE) * (MAX_CU_SIZE);

    for (uint32_t ch = 0; ch < MAX_NUM_COMPONENT; ch++) {
      for (uint32_t buf = 0; buf < 4; buf++) {
        m_yuvExt2[ch][buf] = new Pel[m_yuvExtSize2];
      }
    }
  }

  int shift = bitDepthY + 4;
  for (int i = 32; i < 64; i++) {
    m_auShiftLM[i - 32] = ((1 << shift) + i / 2) / i;
  }
  if (m_piTemp == nullptr) {
    m_piTemp = new Pel[(MAX_CU_SIZE + 1) * (MAX_CU_SIZE + 1)];
  }
  if (m_pMdlmTemp == nullptr) {
    m_pMdlmTemp =
        new Pel[(2 * MAX_CU_SIZE + 1) * (2 * MAX_CU_SIZE + 1)];  // MDLM will use top-above and left-below samples.
  }
}

// ====================================================================================================================
// Public member functions
// ====================================================================================================================

// Function for calculating DC value of the reference samples used in Intra prediction
// NOTE: Bit-Limit - 25-bit source
Pel IntraPrediction::xGetPredValDc(const CPelBuf &pSrc, const Size &dstSize, const int mrlIdx) {
  CHECK(dstSize.width == 0 || dstSize.height == 0, "Empty area provided");

  int idx, sum = 0;
  Pel dcVal;
  const int width = dstSize.width;
  const int height = dstSize.height;
  const auto denom = (width == height) ? (width << 1) : std::max(width, height);
  const auto divShift = getLog2(denom);
  const auto divOffset = (denom >> 1);

  if (width >= height) {
    for (idx = 0; idx < width; idx++) {
      sum += pSrc.at(mrlIdx + 1 + idx, 0);
    }
  }
  if (width <= height) {
    for (idx = 0; idx < height; idx++) {
      sum += pSrc.at(0, mrlIdx + 1 + idx);
    }
  }

  dcVal = (sum + divOffset) >> divShift;
  return dcVal;
}

int IntraPrediction::getWideAngle(int width, int height, int predMode) {
  if (predMode > DC_IDX && predMode <= VDIA_IDX) {
    int modeShift[] = {0, 6, 10, 12, 14, 15};
    int deltaSize = abs(getLog2(width) - getLog2(height));
    if (width > height && predMode < 2 + modeShift[deltaSize]) {
      predMode += (VDIA_IDX - 1);
    } else if (height > width && predMode > VDIA_IDX - modeShift[deltaSize]) {
      predMode -= (VDIA_IDX - 1);
    }
  }
  return predMode;
}

void IntraPrediction::setReferenceArrayLengths(const CompArea &area) {
  // set Top and Left reference samples length
  const int width = area.width;
  const int height = area.height;

  m_leftRefLength = (height << 1);
  m_topRefLength = (width << 1);
}

void IntraPrediction::predIntraAng(const ComponentID compId, PelBuf &piPred, const PredictionUnit &pu,
                                   const bool useFilteredPredSamples) {
  const ComponentID compID = MAP_CHROMA(compId);
  const ChannelType channelType = toChannelType(compID);
  const int iWidth = piPred.width;
  const int iHeight = piPred.height;
  const Size cuSize = Size(pu.blocks[compId].width, pu.blocks[compId].height);
  CHECKD(CU::isMIP(pu, toChannelType(compId)), "We should not get here for MIP.");
  const uint32_t uiDirMode = isLuma(compId) && pu.bdpcmMode()          ? BDPCM_IDX
                             : !isLuma(compId) && pu.bdpcmModeChroma() ? BDPCM_IDX
                                                                       : PU::getFinalIntraMode(pu, channelType);

  CHECKD(iWidth == 2, "Width of 2 is not supported");

  const int multiRefIdx = (compID == COMPONENT_Y) ? pu.multiRefIdx() : 0;
  const bool useISP = pu.ispMode() && isLuma(compID);
  const int srcStride = m_topRefLength + 1 + multiRefIdx;
  const int srcHStride = m_leftRefLength + 1 + multiRefIdx;
  const ClpRng &clpRng(pu.slice->clpRng(compID));
  bool doPDPC = (iWidth >= MIN_TB_SIZEY && iHeight >= MIN_TB_SIZEY) && multiRefIdx == 0;

  const PelBuf &srcBuf = pu.ispMode() && isLuma(compID)
                             ? getISPBuffer(useFilteredPredSamples)
                             : PelBuf(getPredictorPtr(compID, useFilteredPredSamples), srcStride, srcHStride);

  switch (uiDirMode) {
    case (PLANAR_IDX):
      xPredIntraPlanar(srcBuf, piPred, *pu.cs->sps);
      break;
    case (DC_IDX):
      xPredIntraDc(srcBuf, piPred, multiRefIdx);
      break;
    case (BDPCM_IDX):
      xPredIntraBDPCM(srcBuf, piPred, isLuma(compID) ? pu.bdpcmMode() : pu.bdpcmModeChroma(), clpRng);
      break;
    case (2):
    case (DIA_IDX):
    case (VDIA_IDX):
      if (getWideAngle(useISP ? cuSize.width : iWidth, useISP ? cuSize.height : iHeight, uiDirMode) ==
          static_cast<int>(uiDirMode))  // check if uiDirMode is not wide-angle
      {
        xPredIntraAng(srcBuf, piPred, channelType, uiDirMode, clpRng, multiRefIdx, useFilteredPredSamples, doPDPC,
                      useISP, cuSize);
        break;
      }
    default:
      xPredIntraAng(srcBuf, piPred, channelType, uiDirMode, clpRng, multiRefIdx, useFilteredPredSamples, doPDPC, useISP,
                    cuSize);
      break;
  }

  if (doPDPC && (uiDirMode == PLANAR_IDX || uiDirMode == DC_IDX)) {
    if (iWidth > 8)
      IntraPredSampleFilter16(srcBuf.buf, srcBuf.stride, piPred, uiDirMode, clpRng);
    else
      IntraPredSampleFilter8(srcBuf.buf, srcBuf.stride, piPred, uiDirMode, clpRng);
  }
}
void IntraPrediction::predIntraChromaLM(const ComponentID compID, PelBuf &piPred, const PredictionUnit &pu,
                                        const CompArea &chromaArea, int intraDir) {
#if ADAPTIVE_BIT_DEPTH
  if (m_bytePerPixel == 1) {
    xGetLumaRecPixels<Pel8bit>(pu, chromaArea);
    int iLumaStride = 0;
    PelBuf Temp;
    Pel8bit *tempPtr;

    if ((intraDir == MDLM_L_IDX) || (intraDir == MDLM_T_IDX)) {
      iLumaStride = 2 * MAX_CU_SIZE + 1;
      tempPtr = (reinterpret_cast<Pel8bit *>(m_pMdlmTemp)) + iLumaStride + 1;
    } else {
      iLumaStride = MAX_CU_SIZE + 1;
      tempPtr = (reinterpret_cast<Pel8bit *>(m_piTemp)) + iLumaStride + 1;
    }
    Temp = PelBuf(reinterpret_cast<Pel *>(tempPtr), iLumaStride, Size(chromaArea));
    int a, b, iShift;
    xGetLMParameters<Pel8bit>(pu, compID, chromaArea, a, b, iShift);

    ////// final prediction
    piPred.copyFrom2(Temp, 1);
    piPred.linearTransform(a, iShift, b, pu.slice->clpRng(compID));
  } else {
    xGetLumaRecPixels<Pel>(pu, chromaArea);
    int iLumaStride = 0;
    PelBuf Temp;
    if ((intraDir == MDLM_L_IDX) || (intraDir == MDLM_T_IDX)) {
      iLumaStride = 2 * MAX_CU_SIZE + 1;
      Temp = PelBuf(m_pMdlmTemp + iLumaStride + 1, iLumaStride, Size(chromaArea));
    } else {
      iLumaStride = MAX_CU_SIZE + 1;
      Temp = PelBuf(m_piTemp + iLumaStride + 1, iLumaStride, Size(chromaArea));
    }
    int a, b, iShift;
    xGetLMParameters<Pel>(pu, compID, chromaArea, a, b, iShift);

    ////// final prediction
    piPred.copyFrom2(Temp, 2);
    piPred.linearTransform(a, iShift, b, pu.slice->clpRng(compID));
  }
#else
  xGetLumaRecPixels<Pel>(pu, chromaArea);

  int iLumaStride = 0;
  PelBuf Temp;
  if ((intraDir == MDLM_L_IDX) || (intraDir == MDLM_T_IDX)) {
    iLumaStride = 2 * MAX_CU_SIZE + 1;
    Temp = PelBuf(m_pMdlmTemp + iLumaStride + 1, iLumaStride, Size(chromaArea));
  } else {
    iLumaStride = MAX_CU_SIZE + 1;
    Temp = PelBuf(m_piTemp + iLumaStride + 1, iLumaStride, Size(chromaArea));
  }
  int a, b, iShift;
  xGetLMParameters<Pel>(pu, compID, chromaArea, a, b, iShift);

  ////// final prediction
  piPred.copyFrom(Temp);
  piPred.linearTransform(a, iShift, b, pu.slice->clpRng(compID));
#endif
}

void IntraPrediction::xFilterGroup(Pel *pMulDst[], int i, Pel const *const piSrc, int iRecStride, bool bAboveAvaillable,
                                   bool bLeftAvaillable) {
  pMulDst[0][i] = (piSrc[1] + piSrc[iRecStride + 1] + 1) >> 1;

  pMulDst[1][i] = (piSrc[iRecStride] + piSrc[iRecStride + 1] + 1) >> 1;

  pMulDst[3][i] = (piSrc[0] + piSrc[1] + 1) >> 1;

  pMulDst[2][i] = (piSrc[0] + piSrc[1] + piSrc[iRecStride] + piSrc[iRecStride + 1] + 2) >> 2;
}

void IntraPrediction::xPredIntraDc(const CPelBuf &pSrc, PelBuf &pDst, const int mrlIdx) {
#if ADAPTIVE_BIT_DEPTH
  const int width = pDst.width;
  const int height = pDst.height;
  const auto denom = (width == height) ? (width << 1) : std::max(width, height);
  const auto divShift = getLog2(denom);
  const auto divOffset = (denom >> 1);
  int idx, sum = 0;

  if (m_bytePerPixel == 1) {
    const Pel8bit *buf = reinterpret_cast<const Pel8bit *>(pSrc.buf);
    ptrdiff_t srcStride = pSrc.stride;

    if (width >= height) {
      for (idx = 0; idx < width; idx++) {
        sum += buf[mrlIdx + 1 + idx];
      }
    }
    if (width <= height) {
      for (idx = 0; idx < height; idx++) {
        sum += buf[(mrlIdx + 1 + idx) * srcStride];
      }
    }

    Pel8bit dcVal = (sum + divOffset) >> divShift;
    Pel8bit *dstPtr = reinterpret_cast<Pel8bit *>(pDst.buf);
    ptrdiff_t dstStride = pDst.stride;
    if (width == dstStride) {
      std::fill_n(dstPtr, width * height, dcVal);
    } else {
      for (int y = 0; y < height; y++, dstPtr += dstStride) {
        std::fill_n(dstPtr, width, dcVal);
      }
    }
  } else {
    if (width >= height) {
      for (idx = 0; idx < width; idx++) {
        sum += pSrc.at(mrlIdx + 1 + idx, 0);
      }
    }
    if (width <= height) {
      for (idx = 0; idx < height; idx++) {
        sum += pSrc.at(0, mrlIdx + 1 + idx);
      }
    }

    Pel dcVal = (sum + divOffset) >> divShift;
    Pel *dstPtr = pDst.buf;
    ptrdiff_t dstStride = pDst.stride;
    if (width == dstStride) {
      std::fill_n(dstPtr, width * height, dcVal);
    } else {
      for (int y = 0; y < height; y++, dstPtr += dstStride) {
        std::fill_n(dstPtr, width, dcVal);
      }
    }
  }
#else
  const Pel dcval = xGetPredValDc(pSrc, pDst, mrlIdx);
  pDst.fill(dcval);
#endif
}

/** Function for deriving the simplified angular intra predictions.
 *
 * This function derives the prediction samples for the angular mode based on the prediction direction indicated by
 * the prediction mode index. The prediction direction is given by the displacement of the bottom row of the block and
 * the reference row above the block in the case of vertical prediction or displacement of the rightmost column
 * of the block and reference column left from the block in the case of the horizontal prediction. The displacement
 * is signalled at 1/32 pixel accuracy. When projection of the predicted pixel falls inbetween reference samples,
 * the predicted value for the pixel is linearly interpolated from the reference samples. All reference samples are
 * taken from the extended main reference.
 */
// NOTE: Bit-Limit - 25-bit source

void IntraPrediction::xPredIntraAng(const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType,
                                    const uint32_t dirMode, const ClpRng &clpRng, int multiRefIdx,
                                    const bool useFilteredPredSamples, bool &doPDPC, const bool useISP,
                                    const Size cuSize) {
#if ADAPTIVE_BIT_DEPTH
  if (m_bytePerPixel == 1) {
    xPredIntraAngImp<Pel8bit>(pSrc, pDst, channelType, dirMode, clpRng, multiRefIdx, useFilteredPredSamples, doPDPC,
                              useISP, cuSize);
  } else {
    xPredIntraAngImp<Pel>(pSrc, pDst, channelType, dirMode, clpRng, multiRefIdx, useFilteredPredSamples, doPDPC, useISP,
                          cuSize);
  }
#else
  xPredIntraAngImp<Pel>(pSrc, pDst, channelType, dirMode, clpRng, multiRefIdx, useFilteredPredSamples, doPDPC, useISP,
                        cuSize);
#endif
}

template <typename T>
void IntraPrediction::xPredIntraAngImp(const CPelBuf &pSrc, PelBuf &pDst, const ChannelType channelType,
                                       const uint32_t dirMode, const ClpRng &clpRng, int multiRefIdx,
                                       const bool useFilteredPredSamples, bool &doPDPC, const bool useISP,
                                       const Size cuSize) {
  int width = static_cast<int>(pDst.width);
  int height = static_cast<int>(pDst.height);

  CHECK(!(dirMode > DC_IDX && dirMode < NUM_LUMA_MODE), "Invalid intra dir");
  int predMode = useISP ? getWideAngle(cuSize.width, cuSize.height, dirMode) : getWideAngle(width, height, dirMode);
  const bool bIsModeVer = predMode >= DIA_IDX;
  const int intraPredAngleMode = (bIsModeVer) ? predMode - VER_IDX : -(predMode - HOR_IDX);
  const int absAngMode = abs(intraPredAngleMode);
  const int signAng = intraPredAngleMode < 0 ? -1 : 1;

  // Set bitshifts and scale the angle parameter to block size
  static const int angTable[32] = {0,  1,  2,  3,  4,  6,  8,  10, 12, 14,  16,  18,  20,  23,  26,  29,
                                   32, 35, 39, 45, 51, 57, 64, 73, 86, 102, 128, 171, 256, 341, 512, 1024};
  static const int invAngTable[32] = {0,   16384, 8192, 5461, 4096, 2731, 2048, 1638, 1365, 1170, 1024,
                                      910, 819,   712,  630,  565,  512,  468,  420,  364,  321,  287,
                                      256, 224,   191,  161,  128,  96,   64,   48,   32,   16};  // (512 * 32) / Angle
  int invAngle = invAngTable[absAngMode];
  int absAng = angTable[absAngMode];
  int intraPredAngle = signAng * absAng;

  T *refMain;
  T *refSide;

  T refAbove[2 * MAX_CU_SIZE + 3 + 33 * MAX_REF_LINE_IDX];
  T refLeft[2 * MAX_CU_SIZE + 3 + 33 * MAX_REF_LINE_IDX];

  const T *srcPtr = reinterpret_cast<const T *>(pSrc.buf);
  ptrdiff_t srcStride = pSrc.stride;

  // Initialize the Main and Left reference array.
  if (intraPredAngle < 0) {
    for (int x = 0; x <= width + 1 + multiRefIdx; x++) {
      refAbove[x + height] = srcPtr[x];
    }
    for (int y = 0; y <= height + 1 + multiRefIdx; y++) {
      refLeft[y + width] = srcPtr[y * srcStride];
    }
    refMain = bIsModeVer ? refAbove + height : refLeft + width;
    refSide = bIsModeVer ? refLeft + width : refAbove + height;

    // Extend the Main reference to the left.
    int sizeSide = bIsModeVer ? height : width;
    for (int k = -sizeSide; k <= -1; k++) {
      refMain[k] = refSide[std::min((-k * invAngle + 256) >> 9, sizeSide)];
    }
  } else {
    memcpy(refAbove, srcPtr, (m_topRefLength + multiRefIdx + 1) * sizeof(T));
    for (int y = 0; y <= m_leftRefLength + multiRefIdx; y++) {
      refLeft[y] = srcPtr[y * srcStride];
    }

    refMain = bIsModeVer ? refAbove : refLeft;
    refSide = bIsModeVer ? refLeft : refAbove;

    // Extend main reference to right using replication
    const int log2Ratio = getLog2(width) - getLog2(height);
    const int s = std::max<int>(0, bIsModeVer ? log2Ratio : -log2Ratio);
    const int maxIndex = (multiRefIdx << s) + 2;
    const int refLength = bIsModeVer ? m_topRefLength : m_leftRefLength;
    const T val = refMain[refLength + multiRefIdx];
    for (int z = 1; z <= maxIndex; z++) {
      refMain[refLength + multiRefIdx + z] = val;
    }
  }

  // swap width/height if we are doing a horizontal mode:
  Pel tempArray[MAX_TB_SIZEY * MAX_TB_SIZEY];
  const ptrdiff_t dstStride = bIsModeVer ? pDst.stride : MAX_TB_SIZEY;
  Pel *pDstBuf = bIsModeVer ? pDst.buf : tempArray;
  if (!bIsModeVer) {
    std::swap(width, height);
  }

  // compensate for line offset in reference line buffers
  refMain += multiRefIdx;
  refSide += multiRefIdx;

  if (intraPredAngle == 0) {  // pure vertical or pure horizontal
    if (doPDPC) {
      const int scale = ((getLog2(width) - 2 + getLog2(height) - 2 + 2) >> 2);
      CHECK(scale < 0 || scale > 31, "PDPC: scale < 0 || scale > 31");

      IntraHorVerPDPC(pDstBuf, dstStride, reinterpret_cast<Pel *>(refSide), pSrc, width, height, scale,
                      reinterpret_cast<Pel *>(refMain), clpRng);
    } else {
      for (int y = 0; y < height; y++) {
        memcpy(reinterpret_cast<T *>(pDstBuf) + y * dstStride, refMain + 1, width * sizeof(T));
      }
    }
  } else {
    Pel *pDsty = pDstBuf;

    if (!(0 == (absAng & 0x1F))) {
      if (isLuma(channelType)) {
        int deltaPos = intraPredAngle * (1 + multiRefIdx);
        bool interpolationFlag = false, filterFlag = false;
        const int diff = std::min<int>(abs(predMode - HOR_IDX), abs(predMode - VER_IDX));
        const int log2Size = ((getLog2(width) + getLog2(height)) >> 1);
        CHECKD(log2Size >= MAX_INTRA_FILTER_DEPTHS, "Size not supported");
        filterFlag = (diff > m_aucIntraFilter[channelType][log2Size]);

        if (filterFlag) {
          const bool isRefFilter = 0 == (absAng & 0x1F);
          interpolationFlag = !isRefFilter;
        }
        const bool useCubicFilter = useISP ? true : (!interpolationFlag || multiRefIdx > 0);
        const TFilterCoeff *f = (useCubicFilter) ? InterpolationFilter::getChromaFilterTable(0) : g_intraGaussFilter[0];
        if ((width & 7) == 0) {
          IntraPredAngleCore8(pDstBuf, dstStride, reinterpret_cast<Pel *>(refMain), width, height, deltaPos,
                              intraPredAngle, f, useCubicFilter, clpRng);

        } else if ((width & 3) == 0) {
          IntraPredAngleCore4(pDstBuf, dstStride, reinterpret_cast<Pel *>(refMain), width, height, deltaPos,
                              intraPredAngle, f, useCubicFilter, clpRng);
        } else {
          IntraPredAngleCore(pDstBuf, dstStride, reinterpret_cast<Pel *>(refMain), width, height, deltaPos,
                             intraPredAngle, f, useCubicFilter, clpRng);
        }

      } else {
        int deltaPos = intraPredAngle * (1 + multiRefIdx);
        if (width >= 8) {
          IntraPredAngleChroma8(pDstBuf, dstStride, reinterpret_cast<Pel *>(refMain), width, height, deltaPos,
                                intraPredAngle);
        } else if (width == 4) {
          IntraPredAngleChroma4(pDstBuf, dstStride, reinterpret_cast<Pel *>(refMain), width, height, deltaPos,
                                intraPredAngle);
        } else {
          IntraPredAngleChroma(pDstBuf, dstStride, reinterpret_cast<Pel *>(refMain), width, height, deltaPos,
                               intraPredAngle);
        }
      }

    } else {
      T *pDstyTemp = reinterpret_cast<T *>(pDsty);
      for (int y = 0, deltaPos = intraPredAngle * (1 + multiRefIdx); y < height;
           y++, deltaPos += intraPredAngle, pDstyTemp += dstStride) {
        const int deltaInt = deltaPos >> 5;
        // Just copy the integer samples
        memcpy(pDstyTemp, refMain + deltaInt + 1, width * sizeof(T));
      }
    }

    pDsty = pDstBuf;
    int angularScale = 0;
    if (intraPredAngle < 0) {
      doPDPC = false;
    } else if (intraPredAngle > 0) {
      const int sideSize = predMode >= DIA_IDX ? pDst.height : pDst.width;
      const int maxScale = 2;

      angularScale = std::min(maxScale, getLog2(sideSize) - (getLog2(3 * invAngle - 2) - 8));
      doPDPC &= angularScale >= 0;
    }

    if (doPDPC) {
      IntraAnglePDPC(pDsty, dstStride, reinterpret_cast<Pel *>(refSide), width, height, angularScale, invAngle);
    }
  }

  // Flip the block if this is the horizontal mode
  if (!bIsModeVer) {
#if ADAPTIVE_BIT_DEPTH
    pDst.transposedFrom(CPelBuf(pDstBuf, dstStride, width, height), m_bytePerPixel);
#else
    pDst.transposedFrom(CPelBuf(pDstBuf, dstStride, width, height), 2);
#endif
  }
}

void IntraPrediction::xPredIntraBDPCM(const CPelBuf &pSrc, PelBuf &pDst, const uint32_t dirMode, const ClpRng &clpRng) {
  const int wdt = pDst.width;
  const int hgt = pDst.height;

  const ptrdiff_t strideP = pDst.stride;
  const ptrdiff_t strideS = pSrc.stride;

  CHECK(!(dirMode == 1 || dirMode == 2), "Incorrect BDPCM mode parameter.");

#if ADAPTIVE_BIT_DEPTH
  if (m_bytePerPixel == 1) {
    Pel8bit *pred = reinterpret_cast<Pel8bit *>(pDst.buf);
    const Pel8bit *src = reinterpret_cast<const Pel8bit *>(pSrc.buf);
    if (dirMode == 1) {
      for (int y = 0; y < hgt; y++) {
        Pel8bit val = src[(y + 1) * strideS];
        for (int x = 0; x < wdt; x++) {
          pred[x] = val;
        }
        pred += strideP;
      }
    } else {
      for (int y = 0; y < hgt; y++) {
        for (int x = 0; x < wdt; x++) {
          pred[x] = src[x + 1];
        }
        pred += strideP;
      }
    }
  } else {
    Pel *pred = &pDst.buf[0];
    if (dirMode == 1) {
      Pel val;
      for (int y = 0; y < hgt; y++) {
        val = pSrc.buf[(y + 1) * strideS];
        for (int x = 0; x < wdt; x++) {
          pred[x] = val;
        }
        pred += strideP;
      }
    } else {
      for (int y = 0; y < hgt; y++) {
        for (int x = 0; x < wdt; x++) {
          pred[x] = pSrc.buf[x + 1];
        }
        pred += strideP;
      }
    }
  }
#else
  Pel *pred = &pDst.buf[0];
  if (dirMode == 1) {
    Pel val;
    for (int y = 0; y < hgt; y++) {
      val = pSrc.buf[(y + 1) * strideS];
      for (int x = 0; x < wdt; x++) {
        pred[x] = val;
      }
      pred += strideP;
    }
  } else {
    for (int y = 0; y < hgt; y++) {
      for (int x = 0; x < wdt; x++) {
        pred[x] = pSrc.buf[x + 1];
      }
      pred += strideP;
    }
  }
#endif
}

void IntraPrediction::geneWeightedPred(const ComponentID compId, PelBuf &pred, const PredictionUnit &pu, Pel *srcBuf) {
  const int width = pred.width;
  const int height = pred.height;
  const ptrdiff_t srcStride = width;
  const ptrdiff_t dstStride = pred.stride;

  CHECKD(width == 2, "Width of 2 is not supported");

  const CodingUnit &cu = pu;
  const Position posBL = pu.Y().bottomLeft();
  const Position posTR = pu.Y().topRight();

  const CodingUnit *cuLeft = cu.cs->getCURestricted(posBL.offset(-1, 0), cu, CHANNEL_TYPE_LUMA, cu.left);
  const CodingUnit *cuAbove = cu.cs->getCURestricted(posTR.offset(0, -1), cu, CHANNEL_TYPE_LUMA, cu.above);

  const bool isNeigh0Intra = cuLeft && (CU::isIntra(*cuLeft));
  const bool isNeigh1Intra = cuAbove && (CU::isIntra(*cuAbove));

  const int wIntra = 3 - !isNeigh0Intra - !isNeigh1Intra;
  const int wMerge = 3 - !!isNeigh0Intra - !!isNeigh1Intra;

#if ADAPTIVE_BIT_DEPTH
  if (m_bytePerPixel == 1) {
    Pel8bit *dstBuf = reinterpret_cast<Pel8bit *>(pred.buf);
    Pel8bit *srcBufTemp = reinterpret_cast<Pel8bit *>(srcBuf);
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x += 4) {
        dstBuf[y * dstStride + x + 0] =
            (wMerge * dstBuf[y * dstStride + x + 0] + wIntra * srcBufTemp[y * srcStride + x + 0] + 2) >> 2;
        dstBuf[y * dstStride + x + 1] =
            (wMerge * dstBuf[y * dstStride + x + 1] + wIntra * srcBufTemp[y * srcStride + x + 1] + 2) >> 2;
        dstBuf[y * dstStride + x + 2] =
            (wMerge * dstBuf[y * dstStride + x + 2] + wIntra * srcBufTemp[y * srcStride + x + 2] + 2) >> 2;
        dstBuf[y * dstStride + x + 3] =
            (wMerge * dstBuf[y * dstStride + x + 3] + wIntra * srcBufTemp[y * srcStride + x + 3] + 2) >> 2;
      }
    }
  } else {
    Pel *dstBuf = pred.buf;
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x += 4) {
        dstBuf[y * dstStride + x + 0] =
            (wMerge * dstBuf[y * dstStride + x + 0] + wIntra * srcBuf[y * srcStride + x + 0] + 2) >> 2;
        dstBuf[y * dstStride + x + 1] =
            (wMerge * dstBuf[y * dstStride + x + 1] + wIntra * srcBuf[y * srcStride + x + 1] + 2) >> 2;
        dstBuf[y * dstStride + x + 2] =
            (wMerge * dstBuf[y * dstStride + x + 2] + wIntra * srcBuf[y * srcStride + x + 2] + 2) >> 2;
        dstBuf[y * dstStride + x + 3] =
            (wMerge * dstBuf[y * dstStride + x + 3] + wIntra * srcBuf[y * srcStride + x + 3] + 2) >> 2;
      }
    }
  }
#else
  Pel *dstBuf = pred.buf;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x += 4) {
      dstBuf[y * dstStride + x + 0] =
          (wMerge * dstBuf[y * dstStride + x + 0] + wIntra * srcBuf[y * srcStride + x + 0] + 2) >> 2;
      dstBuf[y * dstStride + x + 1] =
          (wMerge * dstBuf[y * dstStride + x + 1] + wIntra * srcBuf[y * srcStride + x + 1] + 2) >> 2;
      dstBuf[y * dstStride + x + 2] =
          (wMerge * dstBuf[y * dstStride + x + 2] + wIntra * srcBuf[y * srcStride + x + 2] + 2) >> 2;
      dstBuf[y * dstStride + x + 3] =
          (wMerge * dstBuf[y * dstStride + x + 3] + wIntra * srcBuf[y * srcStride + x + 3] + 2) >> 2;
    }
  }
#endif
}

void IntraPrediction::switchBuffer(const PredictionUnit &pu, ComponentID compID, PelBuf srcBuff, Pel *dst) {
  Pel *src = srcBuff.bufAt(0, 0);
  int compWidth = compID == COMPONENT_Y ? pu.Y().width : pu.Cb().width;
  int compHeight = compID == COMPONENT_Y ? pu.Y().height : pu.Cb().height;
  for (int i = 0; i < compHeight; i++) {
    memcpy(dst, src, compWidth * sizeof(Pel));
    src += srcBuff.stride;
    dst += compWidth;
  }
}

void IntraPrediction::geneIntrainterPred(const CodingUnit &cu) {
  if (!cu.ciipFlag()) {
    return;
  }
  PROFILER_SCOPE_AND_STAGE_EXT(1, g_timeProfiler, P_INTRAPRED, *cu.cs, compID);

  const PredictionUnit &pu = cu;

  PelUnitBuf predBuf;
  predBuf.bufs.resize(3);

#if JVET_Q0438_MONOCHROME_BUGFIXES
  int maxCompID = 1;
  if (isChromaEnabled(pu.chromaFormat)) {
    maxCompID = MAX_NUM_COMPONENT;
  }
  for (int currCompID = 0; currCompID < maxCompID; currCompID++)
#else
  for (int currCompID = 0; currCompID < 3; currCompID++)
#endif
  {
    if (currCompID > 0 && pu.chromaSize().width <= 2) continue;

    const ComponentID currCompID2 = (ComponentID)currCompID;

    predBuf.bufs[currCompID] = PelBuf(getPredictorPtr2(currCompID2, 0), cu.blocks[currCompID]);
  }

  const bool isUseFilter = IntraPrediction::useFilteredIntraRefSamples(COMPONENT_Y, pu, cu);
  initIntraPatternChType(cu.firstTU, pu.Y(), isUseFilter);
  predIntraAng(COMPONENT_Y, predBuf.Y(), pu, isUseFilter);

#if JVET_Q0438_MONOCHROME_BUGFIXES
  if (isChromaEnabled(pu.chromaFormat) && pu.chromaSize().width > 2)
#else
  if (pu.chromaSize().width > 2)
#endif
  {
    initIntraPatternChType(cu.firstTU, pu.Cb(), false);
    predIntraAng(COMPONENT_Cb, predBuf.Cb(), pu, false);

    initIntraPatternChType(cu.firstTU, pu.Cr(), false);
    predIntraAng(COMPONENT_Cr, predBuf.Cr(), pu, false);
  }
}

inline int isAboveAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posLT,
                            const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, int &neighborSize);
inline int isLeftAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posLT,
                           const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, int &neighborSize);

inline int isAboveRightAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posRT,
                                 const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, int &neighborSize);
inline int isBelowLeftAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posLB,
                                const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, int &neighborSize);

void IntraPrediction::initIntraPatternChType(const TransformUnit &tu, const CompArea &area,
                                             const bool bFilterRefSamples) {
  CHECK(area.width == 2, "Width of 2 is not supported");
  const CodingStructure &cs = *tu.cu->cs;

  Pel *refBufUnfiltered = m_piYuvExt[area.compID][PRED_BUF_UNFILTERED];
  Pel *refBufFiltered = m_piYuvExt[area.compID][PRED_BUF_FILTERED];

  setReferenceArrayLengths(area);

#if ADAPTIVE_BIT_DEPTH
  if (m_bytePerPixel == 1) {
    // ----- Step 1: unfiltered reference samples -----
    xFillReferenceSamples<Pel8bit>(cs.picture->getRecoBuf(area, false, m_bytePerPixel), refBufUnfiltered, area, tu);
    // ----- Step 2: filtered reference samples -----
    if (bFilterRefSamples) {
      xFilterReferenceSamples<Pel8bit>(refBufUnfiltered, refBufFiltered, area, *cs.sps, tu.cu->multiRefIdx());
    }
  } else {
    // ----- Step 1: unfiltered reference samples -----
    xFillReferenceSamples<Pel>(cs.picture->getRecoBuf(area, false, m_bytePerPixel), refBufUnfiltered, area, tu);
    // ----- Step 2: filtered reference samples -----
    if (bFilterRefSamples) {
      xFilterReferenceSamples<Pel>(refBufUnfiltered, refBufFiltered, area, *cs.sps, tu.cu->multiRefIdx());
    }
  }
#else
  // ----- Step 1: unfiltered reference samples -----
  xFillReferenceSamples<Pel>(cs.picture->getRecoBuf(area), refBufUnfiltered, area, tu);
  // ----- Step 2: filtered reference samples -----
  if (bFilterRefSamples) {
    xFilterReferenceSamples<Pel>(refBufUnfiltered, refBufFiltered, area, *cs.sps, tu.cu->multiRefIdx());
  }
#endif
}

void IntraPrediction::initIntraPatternChTypeISP(const CodingUnit &cu, const CompArea &area, PelBuf &recBuf) {
  const CodingStructure &cs = *cu.cs;

  const Position &posLT = area.pos();
  bool isLeftAvail = nullptr != cs.getCURestricted(posLT.offset(-1, 0), cu, CH_L, posLT.x == cu.lx() ? cu.left : &cu);
  bool isAboveAvail = nullptr != cs.getCURestricted(posLT.offset(0, -1), cu, CH_L, posLT.y == cu.ly() ? cu.left : &cu);

  // ----- Step 1: unfiltered reference samples -----
  if (cu.blocks[area.compID].x == area.x && cu.blocks[area.compID].y == area.y) {
    Pel *refBufUnfiltered = m_piYuvExt[area.compID][PRED_BUF_UNFILTERED];
    if (cu.ispMode() == HOR_INTRA_SUBPARTITIONS) {
      m_leftRefLength = cu.Y().height << 1;
      m_topRefLength = cu.Y().width + area.width;
    } else  // if (cu.ispMode() == VER_INTRA_SUBPARTITIONS)
    {
      m_leftRefLength = cu.Y().height + area.height;
      m_topRefLength = cu.Y().width << 1;
    }

    const int srcStride = m_topRefLength + 1;
    const int srcHStride = m_leftRefLength + 1;

    m_pelBufISP[0] = m_pelBufISPBase[0] = PelBuf(m_piYuvExt[area.compID][PRED_BUF_UNFILTERED], srcStride, srcHStride);
    m_pelBufISP[1] = m_pelBufISPBase[1] = PelBuf(m_piYuvExt[area.compID][PRED_BUF_FILTERED], srcStride, srcHStride);

#if ADAPTIVE_BIT_DEPTH
    if (m_bytePerPixel == 1) {
      xFillReferenceSamples<Pel8bit>(cs.picture->getRecoBuf(cu.Y(), false, m_bytePerPixel), refBufUnfiltered, cu.Y(),
                                     isLuma(area.compID) ? cu.firstTU : *cu.lastTU);
    } else {
      xFillReferenceSamples<Pel>(cs.picture->getRecoBuf(cu.Y(), false, m_bytePerPixel), refBufUnfiltered, cu.Y(),
                                 isLuma(area.compID) ? cu.firstTU : *cu.lastTU);
    }
#else
    xFillReferenceSamples<Pel>(cs.picture->getRecoBuf(cu.Y()), refBufUnfiltered, cu.Y(),
                               isLuma(area.compID) ? cu.firstTU : *cu.lastTU);
#endif

    // After having retrieved all the CU reference samples, the number of reference samples is now adjusted for the
    // current subpartition
    m_topRefLength = cu.blocks[area.compID].width + area.width;
    m_leftRefLength = cu.blocks[area.compID].height + area.height;
  } else {
    // Now we only need to fetch the newly available reconstructed samples from the previously coded TU
#if ADAPTIVE_BIT_DEPTH
    if (m_bytePerPixel == 1) {
      xFillReferenceSamplesSameCu<Pel8bit>(area, cu, recBuf, isLeftAvail, isAboveAvail);
    } else {
      xFillReferenceSamplesSameCu<Pel>(area, cu, recBuf, isLeftAvail, isAboveAvail);
    }
#else
    xFillReferenceSamplesSameCu<Pel>(area, cu, recBuf, isLeftAvail, isAboveAvail);
#endif
  }
}

template <typename T>
void IntraPrediction::xFillReferenceSamples(const CPelBuf &recoBuf, Pel *_refBufUnfiltered, const CompArea &area,
                                            const TransformUnit &tu) const {
  const ChannelType chType = toChannelType(area.compID);
  const CodingUnit &cu = *tu.cu;
  const CodingStructure &cs = *cu.cs;
  const SPS &sps = *cs.sps;
  const PreCalcValues &pcv = *cs.pcv;

  const int multiRefIdx = (area.compID == COMPONENT_Y) ? cu.multiRefIdx() : 0;

  const int tuWidth = area.width;
  const int tuHeight = area.height;
  const int predSize = m_topRefLength;
  const int predHSize = m_leftRefLength;
  const int predStride = predSize + 1 + multiRefIdx;

  const int unitWidth = tuWidth <= 2 && cu.ispMode() && isLuma(area.compID)
                            ? tuWidth
                            : pcv.minCUWidth >> getComponentScaleX(area.compID, sps.getChromaFormatIdc());
  const int unitHeight = tuHeight <= 2 && cu.ispMode() && isLuma(area.compID)
                             ? tuHeight
                             : pcv.minCUHeight >> getComponentScaleY(area.compID, sps.getChromaFormatIdc());

  const int totalAboveUnits = (predSize + (unitWidth - 1)) / unitWidth;
  const int totalLeftUnits = (predHSize + (unitHeight - 1)) / unitHeight;
  const int totalUnits = totalAboveUnits + totalLeftUnits + 1;  //+1 for top-left
  const int numAboveUnits = std::max<int>(tuWidth / unitWidth, 1);
  const int numLeftUnits = std::max<int>(tuHeight / unitHeight, 1);
  const int numAboveRightUnits = totalAboveUnits - numAboveUnits;
  const int numLeftBelowUnits = totalLeftUnits - numLeftUnits;

  CHECKD(numAboveUnits <= 0 || numLeftUnits <= 0 || numAboveRightUnits <= 0 || numLeftBelowUnits <= 0,
         "Size not supported");

  // ----- Step 1: analyze neighborhood -----
  const Position posLT = area.pos();

  int neighborSize[3] = {0};

  const CodingUnit *aboveLeftCu =
      cu.cs->getCURestricted(posLT.offset(-1, -1), cu, chType, cu.left ? cu.left : cu.above);
  neighborSize[0] = !!aboveLeftCu;
  int numIntraNeighbor = aboveLeftCu ? 1 : 0;
  numIntraNeighbor +=
      isAboveAvailable(tu, chType, posLT, numAboveUnits + numAboveRightUnits, unitWidth, neighborSize[1]);
  numIntraNeighbor += isLeftAvailable(tu, chType, posLT, numLeftUnits + numLeftBelowUnits, unitHeight, neighborSize[2]);

  // ----- Step 2: fill reference samples (depending on neighborhood) -----
  CHECK((predHSize + 1) * predStride > m_iYuvExtSize, "Reference sample area not supported");

  T *refBufUnfiltered = reinterpret_cast<T *>(_refBufUnfiltered);
  const T *srcBuf = reinterpret_cast<const T *>(recoBuf.buf);
  const ptrdiff_t srcStride = recoBuf.stride;
  T *ptrDst = refBufUnfiltered;
  const T *ptrSrc;
  const T valueDC = 1 << (sps.getBitDepth(chType) - 1);

  if (numIntraNeighbor == 0) {
    // Fill border with DC value
    std::fill_n(ptrDst, (predSize + multiRefIdx + 1), valueDC);
    for (int i = 1; i <= predHSize + multiRefIdx; i++) {
      ptrDst[i * predStride] = valueDC;
    }
  } else if (numIntraNeighbor == totalUnits) {
    // Fill top-left border and top and top right with rec. samples
    ptrSrc = srcBuf - (1 + multiRefIdx) * srcStride - (1 + multiRefIdx);
    memcpy(ptrDst, ptrSrc, sizeof(T) * (predSize + multiRefIdx + 1));
    ptrSrc = srcBuf - multiRefIdx * srcStride - (1 + multiRefIdx);
    for (int i = 1; i <= predHSize + multiRefIdx; i++) {
      ptrDst[i * predStride] = *(ptrSrc);
      ptrSrc += srcStride;
    }
  } else {  // reference samples are partially available
    // Fill top-left sample(s) if available
    if (neighborSize[2] > 0) {  // left is available
      // Fill left & below-left samples if available (downwards)
      ptrSrc = srcBuf - (1 + multiRefIdx);
      ptrDst = refBufUnfiltered + (1 + multiRefIdx) * predStride;
      int tmpSize = neighborSize[2] * unitHeight;
      tmpSize = std::min(tmpSize, predHSize);
      for (int i = 0; i < tmpSize; i++) {
        ptrDst[i * predStride] = ptrSrc[i * srcStride];
      }

      // pad
      T tmpPixel = ptrDst[(tmpSize - 1) * predStride];
      for (int i = tmpSize; i < predHSize; i++) {
        ptrDst[i * predStride] = tmpPixel;
      }

      // Fill top-left sample(s) if available
      if (neighborSize[0]) {
        ptrSrc = srcBuf - (1 + multiRefIdx) * srcStride - (1 + multiRefIdx);
        ptrDst = refBufUnfiltered;
        memcpy(ptrDst, ptrSrc, sizeof(T) * (multiRefIdx + 1));
        for (int i = 1; i <= multiRefIdx; i++) {
          ptrDst[i * predStride] = ptrSrc[i * srcStride];
        }
      } else {                                // pad
        ptrSrc = srcBuf - (1 + multiRefIdx);  // left pixel
        ptrDst = refBufUnfiltered;
        tmpPixel = ptrSrc[0];
        ptrDst[0] = tmpPixel;
        for (int i = 1; i <= multiRefIdx; i++) {
          ptrDst[i] = tmpPixel;
          ptrDst[i * predStride] = tmpPixel;
        }
      }

      // Fill above & above-right samples if available (left-to-right)
      if (neighborSize[1]) {
        ptrSrc = srcBuf - srcStride * (1 + multiRefIdx);
        ptrDst = refBufUnfiltered + 1 + multiRefIdx;
        tmpSize = neighborSize[1] * unitWidth;
        tmpSize = std::min(tmpSize, predSize);
        memcpy(ptrDst, ptrSrc, tmpSize * sizeof(T));
        // pad
        T tmpPixel = ptrDst[tmpSize - 1];
        for (int i = tmpSize; i < predSize; i++) {
          ptrDst[i] = tmpPixel;
        }
      } else {  // all not available, pad
        ptrSrc = srcBuf - srcStride * (1 + multiRefIdx);
        ptrDst = refBufUnfiltered + 1 + multiRefIdx;
        T tmpPixel = ptrDst[-1];
        std::fill_n(ptrDst, predSize, tmpPixel);
      }
    } else {  // left is not available, top must be available
      // Fill above & above-right samples (left-to-right)
      ptrSrc = srcBuf - srcStride * (1 + multiRefIdx);
      ptrDst = refBufUnfiltered + 1 + multiRefIdx;
      int tmpSize = neighborSize[1] * unitWidth;
      tmpSize = std::min(tmpSize, predSize);
      memcpy(ptrDst, ptrSrc, tmpSize * sizeof(T));
      // pad
      T tmpPixel = ptrDst[tmpSize - 1];
      for (int i = tmpSize; i < predSize; i++) {
        ptrDst[i] = tmpPixel;
      }

      tmpPixel = ptrSrc[0];
      // pad top-left sample(s)
      ptrDst = refBufUnfiltered;
      ptrDst[0] = tmpPixel;
      for (int i = 1; i <= multiRefIdx; i++) {
        ptrDst[i] = tmpPixel;
        ptrDst[i * predStride] = tmpPixel;
      }

      // pad left sample(s)
      ptrDst = refBufUnfiltered + (1 + multiRefIdx) * predStride;
      for (int i = 0; i < predHSize; i++) {
        ptrDst[i * predStride] = tmpPixel;
      }
    }
  }
}

template <typename T>
void IntraPrediction::xFillReferenceSamplesSameCu(const CompArea &area, const CodingUnit &cu, PelBuf &recBuf,
                                                  bool isLeftAvail, bool isAboveAvail) {
  // Now we only need to fetch the newly available reconstructed samples from the previously coded TU
  Position tuPos = area;
  tuPos.relativeTo(cu.Y());

  T *buf0 = (reinterpret_cast<T *>(m_pelBufISPBase[0].buf)) + tuPos.x + tuPos.y * m_pelBufISPBase[0].stride;
  T *buf1 = (reinterpret_cast<T *>(m_pelBufISPBase[1].buf)) + tuPos.x + tuPos.y * m_pelBufISPBase[1].stride;
  m_pelBufISP[0] = AreaBuf<Pel>(reinterpret_cast<Pel *>(buf0), m_pelBufISPBase[0].stride, area.size());
  m_pelBufISP[1] = AreaBuf<Pel>(reinterpret_cast<Pel *>(buf1), m_pelBufISPBase[1].stride, area.size());

  PelBuf &dstBuf = m_pelBufISP[0];

  m_topRefLength = cu.blocks[area.compID].width + area.width;
  m_leftRefLength = cu.blocks[area.compID].height + area.height;

  T *recPtr = reinterpret_cast<T *>(recBuf.buf);
  T *dstPtr = reinterpret_cast<T *>(dstBuf.buf);
  ptrdiff_t recStride = recBuf.stride;
  ptrdiff_t dstStride = dstBuf.stride;

  const int predSizeHor = m_topRefLength;
  const int predSizeVer = m_leftRefLength;
  if (cu.ispMode() == HOR_INTRA_SUBPARTITIONS) {
    T *src = recPtr - recStride;
    T *dst = dstPtr + 1;
    for (int i = 0; i < area.width; i++) {
      dst[i] = src[i];
    }
    T sample = src[area.width - 1];
    dst += area.width;
    for (int i = 0; i < predSizeHor - area.width; i++) {
      dst[i] = sample;
    }
    if (!isLeftAvail) {  // if left is not avaible, then it is necessary to fetch these samples for each subpartition
      T *dst = dstPtr;
      T sample = src[0];
      for (int i = 0; i < predSizeVer + 1; i++) {
        *dst = sample;
        dst += dstStride;
      }
    }
  } else {
    T *src = recPtr - 1;
    T *dst = dstPtr + dstStride;
    for (int i = 0; i < area.height; i++) {
      *dst = *src;
      src += recStride;
      dst += dstStride;
    }
    T sample = src[-recStride];
    for (int i = 0; i < predSizeVer - area.height; i++) {
      *dst = sample;
      dst += dstStride;
    }

    if (!isAboveAvail) {  // if above is not avaible, then it is necessary to fetch these samples for each subpartition
      T *dst = dstPtr;
      T sample = recPtr[-1];
      for (int i = 0; i < predSizeHor + 1; i++) {
        dst[i] = sample;
      }
    }
  }
}

template <typename T>
void IntraPrediction::xFilterReferenceSamples(const Pel *_refBufUnfiltered, Pel *_refBufFiltered, const CompArea &area,
                                              const SPS &sps, int multiRefIdx, ptrdiff_t stride) const {
  const T *refBufUnfiltered = reinterpret_cast<const T *>(_refBufUnfiltered);
  T *refBufFiltered = reinterpret_cast<T *>(_refBufFiltered);

  if (area.compID != COMPONENT_Y) {
    multiRefIdx = 0;
  }
  const int predSize = m_topRefLength + multiRefIdx;
  const int predHSize = m_leftRefLength + multiRefIdx;
  const ptrdiff_t predStride = stride == 0 ? predSize + 1 : stride;

  // Regular reference sample filter
  const T *piSrcPtr = refBufUnfiltered + (predStride * predHSize);  // bottom left
  T *piDestPtr = refBufFiltered + (predStride * predHSize);         // bottom left

  // bottom left (not filtered)
  *piDestPtr = *piSrcPtr;
  piDestPtr -= predStride;
  piSrcPtr -= predStride;
  // left column (bottom to top)
  for (int i = 1; i < predHSize; i++, piDestPtr -= predStride, piSrcPtr -= predStride) {
    *piDestPtr = (piSrcPtr[predStride] + 2 * piSrcPtr[0] + piSrcPtr[-predStride] + 2) >> 2;
  }
  // top-left
  *piDestPtr = (piSrcPtr[predStride] + 2 * piSrcPtr[0] + piSrcPtr[1] + 2) >> 2;
  piDestPtr++;
  piSrcPtr++;
  // top row (left-to-right)
  for (uint32_t i = 1; i < predSize; i++, piDestPtr++, piSrcPtr++) {
    *piDestPtr = (piSrcPtr[1] + 2 * piSrcPtr[0] + piSrcPtr[-1] + 2) >> 2;
  }
  // top right (not filtered)
  *piDestPtr = *piSrcPtr;
}

bool IntraPrediction::getUseFilterRef(const int predMode, const int dirMode) {
  static const int angTable[32] = {0,  1,  2,  3,  4,  6,  8,  10, 12, 14,  16,  18,  20,  23,  26,  29,
                                   32, 35, 39, 45, 51, 57, 64, 73, 86, 102, 128, 171, 256, 341, 512, 1024};

  const int intraPredAngleMode = (predMode >= DIA_IDX) ? predMode - VER_IDX : -(predMode - HOR_IDX);

  const int absAngMode = abs(intraPredAngleMode);
  const int absAng = angTable[absAngMode];

  return 0 == (absAng & 0x1F);
}

bool IntraPrediction::useFilteredIntraRefSamples(const ComponentID &compID, const PredictionUnit &pu,
                                                 const UnitArea &tuArea) {
  // const SPS         &sps    = *pu.cs->sps;
  const ChannelType chType = toChannelType(compID);

  // high level conditions
  // if( sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag() )  { return false; }
  // if( !isLuma( chType ) )                                           { return false; }
  // if( pu.ispMode() && isLuma(compID) )                              { return false; }
  // if( CU::isMIP( pu, chType ) )                                     { return false; }
  if (pu.multiRefIdx()) {
    return false;
  }
  if (pu.bdpcmMode()) {
    return false;
  }

  // pred. mode related conditions
  const int dirMode = PU::getFinalIntraMode(pu, chType);
  int predMode = getWideAngle(tuArea.blocks[compID].width, tuArea.blocks[compID].height, dirMode);
  if (dirMode == DC_IDX) {
    return false;
  }
  if (dirMode == PLANAR_IDX) {
    return tuArea.blocks[compID].area() > 32 ? true : false;
  }

  bool filterFlag = false;
  {
    const int diff = std::min<int>(abs(predMode - HOR_IDX), abs(predMode - VER_IDX));
    int log2Size = ((getLog2(tuArea.blocks[compID].width) + getLog2(tuArea.blocks[compID].height)) >> 1);
    CHECKD(log2Size >= MAX_INTRA_FILTER_DEPTHS, "Size not supported");
    filterFlag = (diff > m_aucIntraFilter[chType][log2Size]);
  }

  if (filterFlag) {
    const bool isRefFilter = getUseFilterRef(predMode, dirMode);
    CHECKD(tuArea.blocks[compID].width * tuArea.blocks[compID].height <= 32,
           "DCT-IF interpolation filter is always used for 4x4, 4x8, and 8x4 luma CB");
    return isRefFilter;
  } else {
    return false;
  }
}

static inline TransformUnit const *getTU(const CodingUnit &cu, const Position &pos, const ChannelType chType) {
  const TransformUnit *ptu = &cu.firstTU;

  if (!ptu->next) return ptu;

  while (!(ptu->blocks[chType].x + ptu->blocks[chType].width > pos.x &&
           ptu->blocks[chType].y + ptu->blocks[chType].height > pos.y)) {
    ptu = ptu->next;
  }

  return ptu;
}

int isAboveAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posLT,
                     const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, int &neighborSize) {
  const CodingUnit &cu = *tu.cu;
  const CodingStructure &cs = *cu.cs;
  neighborSize = 0;
  int maxDx = uiNumUnitsInPU * unitWidth;
  int rightXAbove = -1;
  Position refPos = posLT.offset(0, -1);
  const TransformUnit *pcTUAbove = nullptr;
  int currTUIdx = tu.idx;
  for (uint32_t dx = 0; dx < maxDx; dx += unitWidth, refPos.x += unitWidth) {
    if (!pcTUAbove || refPos.x >= rightXAbove) {
      const CodingUnit *cuAbove = cs.getCURestricted(refPos, cu, chType, pcTUAbove ? nullptr : cu.above);
      if (!cuAbove) break;
      pcTUAbove = getTU(*cuAbove, refPos, chType);
      if (pcTUAbove->idx >= currTUIdx) break;
      rightXAbove = pcTUAbove->blocks[chType].x + pcTUAbove->blocks[chType].width;
    }
    neighborSize++;
  }

  return neighborSize;
}

int isLeftAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posLT,
                    const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, int &neighborSize) {
  const CodingUnit &cu = *tu.cu;
  const CodingStructure &cs = *cu.cs;
  neighborSize = 0;
  int maxDy = uiNumUnitsInPU * unitHeight;
  int bottomYLeft = -1;
  Position refPos = posLT.offset(-1, 0);
  const TransformUnit *pcTULeft = nullptr;
  int currTUIdx = tu.idx;
  for (uint32_t dy = 0; dy < maxDy; dy += unitHeight, refPos.y += unitHeight) {
    if (!pcTULeft || refPos.y >= bottomYLeft) {
      const CodingUnit *cuLeft = cs.getCURestricted(refPos, cu, chType, pcTULeft ? nullptr : cu.left);

      if (!cuLeft) break;
      pcTULeft = getTU(*cuLeft, refPos, chType);
      if (pcTULeft->idx >= currTUIdx) break;
      bottomYLeft = pcTULeft->blocks[chType].y + pcTULeft->blocks[chType].height;
    }
    neighborSize++;
  }
  return neighborSize;
}

int isAboveRightAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posRT,
                          const uint32_t uiNumUnitsInPU, const uint32_t unitWidth, int &neighborSize) {
  const CodingUnit &cu = *tu.cu;
  const CodingStructure &cs = *cu.cs;

  neighborSize = 0;
  int maxDx = uiNumUnitsInPU * unitWidth;
  int rightXAbove = -1;
  Position refPos = posRT.offset(unitWidth, -1);
  const TransformUnit *pcTUAbove = nullptr;
  int currTUIdx = tu.idx;
  for (uint32_t dx = 0; dx < maxDx; dx += unitWidth, refPos.x += unitWidth) {
    if (!pcTUAbove || refPos.x >= rightXAbove) {
      const CodingUnit *cuAbove = cs.getCURestricted(refPos, cu, chType, pcTUAbove ? nullptr : cu.above);

      if (!cuAbove) break;

      pcTUAbove = getTU(*cuAbove, refPos, chType);
      if (pcTUAbove->idx >= currTUIdx) break;
      rightXAbove = pcTUAbove->blocks[chType].x + pcTUAbove->blocks[chType].width;
    }
    neighborSize++;
  }
  return neighborSize;
}

int isBelowLeftAvailable(const TransformUnit &tu, const ChannelType &chType, const Position &posLB,
                         const uint32_t uiNumUnitsInPU, const uint32_t unitHeight, int &neighborSize) {
  const CodingUnit &cu = *tu.cu;
  const CodingStructure &cs = *cu.cs;
  neighborSize = 0;
  int maxDy = uiNumUnitsInPU * unitHeight;
  int bottomYLeft = -1;
  Position refPos = posLB.offset(-1, unitHeight);
  const TransformUnit *pcTULeft = nullptr;
  int currTUIdx = tu.idx;
  for (uint32_t dy = 0; dy < maxDy; dy += unitHeight, refPos.y += unitHeight) {
    if (!pcTULeft || refPos.y >= bottomYLeft) {
      const CodingUnit *cuLeft = cs.getCURestricted(refPos, cu, chType, pcTULeft ? nullptr : cu.left);

      if (!cuLeft) break;

      pcTULeft = getTU(*cuLeft, refPos, chType);
      if (pcTULeft->idx >= currTUIdx) break;
      bottomYLeft = pcTULeft->blocks[chType].y + pcTULeft->blocks[chType].height;
    }
    neighborSize++;
  }

  return neighborSize;
}
// LumaRecPixels
template <typename T>
void IntraPrediction::xGetLumaRecPixels(const PredictionUnit &pu, CompArea chromaArea) {
  int iDstStride = 0;
  T *pDst0 = 0;
  int curChromaMode = pu.intraDir[1];
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX)) {
    iDstStride = 2 * MAX_CU_SIZE + 1;
    pDst0 = (reinterpret_cast<T *>(m_pMdlmTemp)) + iDstStride + 1;
  } else {
    iDstStride = MAX_CU_SIZE + 1;
    pDst0 = (reinterpret_cast<T *>(m_piTemp)) + iDstStride + 1;  // MMLM_SAMPLE_NEIGHBOR_LINES;
  }
  // assert 420 chroma subsampling
  CompArea lumaArea = CompArea(COMPONENT_Y, chromaArea.lumaPos(pu.chromaFormat),
                               recalcSize(pu.chromaFormat, CHANNEL_TYPE_CHROMA, CHANNEL_TYPE_LUMA,
                                          chromaArea.size()));  // needed for correct pos/size (4x4 Tus)

  CHECKD(lumaArea.width == chromaArea.width && CHROMA_444 != pu.chromaFormat, "");
  CHECKD(lumaArea.height == chromaArea.height && CHROMA_444 != pu.chromaFormat && CHROMA_422 != pu.chromaFormat, "");

  const SizeType uiCWidth = chromaArea.width;
  const SizeType uiCHeight = chromaArea.height;
#if ADAPTIVE_BIT_DEPTH
  CPelBuf src = pu.cs->picture->getRecoBuf(lumaArea, false, m_bytePerPixel);
#else
  CPelBuf src = pu.cs->picture->getRecoBuf(lumaArea);
#endif
  T const *pRecSrc0 = (T const *)src.buf;
  ptrdiff_t iRecStride = src.stride;
  int logSubWidthC = getChannelTypeScaleX(CHANNEL_TYPE_CHROMA, pu.chromaFormat);
  int logSubHeightC = getChannelTypeScaleY(CHANNEL_TYPE_CHROMA, pu.chromaFormat);

  ptrdiff_t iRecStride2 = iRecStride << logSubHeightC;  // TODO: really Height here? not Width?
  const int mult = 1 << logSubWidthC;
  const CodingUnit &lumaCU = isChroma(pu.chType()) ? *pu.cs->getCU(lumaArea.pos(), CH_L) : pu;
  const CodingUnit &cu = pu;

  const CompArea &area = isChroma(pu.chType()) ? chromaArea : lumaArea;

  const uint32_t uiTuWidth = area.width;
  const uint32_t uiTuHeight = area.height;

  int iBaseUnitSize = (1 << MIN_CU_LOG2);

  const int iUnitWidth = iBaseUnitSize >> getComponentScaleX(area.compID, pu.chromaFormat);
  const int iUnitHeight = iBaseUnitSize >> getComponentScaleY(area.compID, pu.chromaFormat);
  const int iTUWidthInUnits = uiTuWidth / iUnitWidth;
  const int iTUHeightInUnits = uiTuHeight / iUnitHeight;
  const int iAboveUnits = iTUWidthInUnits;
  const int iLeftUnits = iTUHeightInUnits;
  const int chromaUnitWidth = iBaseUnitSize >> getComponentScaleX(COMPONENT_Cb, pu.chromaFormat);
  const int chromaUnitHeight = iBaseUnitSize >> getComponentScaleY(COMPONENT_Cb, pu.chromaFormat);
  const int topTemplateSampNum = 2 * uiCWidth;  // for MDLM, the number of template samples is 2W or 2H.
  const int leftTemplateSampNum = 2 * uiCHeight;
  CHECKD(!(m_topRefLength >= topTemplateSampNum), "Error!");
  CHECKD(!(m_leftRefLength >= leftTemplateSampNum), "Error!");
  int totalAboveUnits =
      (curChromaMode == MDLM_T_IDX) ? (topTemplateSampNum + (chromaUnitWidth - 1)) / chromaUnitWidth : iAboveUnits;
  int totalLeftUnits =
      (curChromaMode == MDLM_L_IDX) ? (leftTemplateSampNum + (chromaUnitHeight - 1)) / chromaUnitHeight : iLeftUnits;
  int neighborSize;

  const CodingUnit &chromaCU = isChroma(pu.chType()) ? cu : lumaCU;
  const TransformUnit &chromaTU = *getTU(chromaCU, chromaArea.pos(), CH_C);

  int availlableLeftUnit =
      isLeftAvailable(chromaTU, toChannelType(area.compID), area.pos(), totalLeftUnits, iUnitHeight, neighborSize);

  const bool bLeftAvaillable = availlableLeftUnit >= iTUHeightInUnits;

  int availlableAboveUnit =
      isAboveAvailable(chromaTU, toChannelType(area.compID), area.pos(), totalAboveUnits, iUnitWidth, neighborSize);

  const bool bAboveAvaillable = availlableAboveUnit >= iTUWidthInUnits;

  T *pDst = nullptr;
  T const *piSrc = nullptr;

  bool isFirstRowOfCtu = (lumaArea.y & ((pu.cs->sps)->getCTUSize() - 1)) == 0;
  const ptrdiff_t strOffset = (CHROMA_444 == pu.chromaFormat) ? 0 : iRecStride;

  int c0_3tap = 2, c1_3tap = 1, c2_3tap = 1, offset_3tap = 2, shift_3tap = 2;                            // sum = 4
  int c0_5tap = 1, c1_5tap = 4, c2_5tap = 1, c3_5tap = 1, c4_5tap = 1, offset_5tap = 4, shift_5tap = 3;  // sum = 8
  int c0_6tap = 2, c1_6tap = 1, c2_6tap = 1, c3_6tap = 2, c4_6tap = 1, c5_6tap = 1, offset_6tap = 4,
      shift_6tap = 3;  // sum = 8

  switch (pu.chromaFormat) {
    case CHROMA_422:  // overwrite filter coefficient values for 422
      c0_3tap = 2, c1_3tap = 1, c2_3tap = 1, offset_3tap = 2, shift_3tap = 2;                            // sum = 4
      c0_5tap = 0, c1_5tap = 1, c2_5tap = 0, c3_5tap = 0, c4_5tap = 0, offset_5tap = 0, shift_5tap = 0;  // sum = 1
      c0_6tap = 2, c1_6tap = 1, c2_6tap = 1, c3_6tap = 0, c4_6tap = 0, c5_6tap = 0, offset_6tap = 2,
      shift_6tap = 2;  // sum = 4
      break;

    case CHROMA_444:  // overwrite filter coefficient values for 444
      c0_3tap = 1, c1_3tap = 0, c2_3tap = 0, offset_3tap = 0, shift_3tap = 0;                            // sum = 1
      c0_5tap = 0, c1_5tap = 1, c2_5tap = 0, c3_5tap = 0, c4_5tap = 0, offset_5tap = 0, shift_5tap = 0;  // sum = 1
      c0_6tap = 1, c1_6tap = 0, c2_6tap = 0, c3_6tap = 0, c4_6tap = 0, c5_6tap = 0, offset_6tap = 0,
      shift_6tap = 0;  // sum = 1
      break;

    default:
      break;
  }

  if (bAboveAvaillable) {
    pDst = pDst0 - iDstStride;
    int avaiAboveSizes = availlableAboveUnit * chromaUnitWidth;
    for (int i = 0; i < avaiAboveSizes; i++) {
      if (isFirstRowOfCtu) {
        piSrc = pRecSrc0 - iRecStride;

        if ((i == 0 && !bLeftAvaillable) || (i == avaiAboveSizes - 1 + logSubWidthC)) {
          pDst[i] =
              (piSrc[mult * i] * c0_3tap + piSrc[mult * i] * c1_3tap + piSrc[mult * i + 1] * c2_3tap + offset_3tap) >>
              shift_3tap;
        } else {
          pDst[i] = (piSrc[mult * i] * c0_3tap + piSrc[mult * i - 1] * c1_3tap + piSrc[mult * i + 1] * c2_3tap +
                     offset_3tap) >>
                    shift_3tap;
        }
      } else if (pu.cs->sps->getCclmCollocatedChromaFlag()) {
        piSrc = pRecSrc0 - iRecStride2;

        if ((i == 0 && !bLeftAvaillable) || (i == avaiAboveSizes - 1 + logSubWidthC)) {
          pDst[i] = (piSrc[mult * i - strOffset] * c0_5tap + piSrc[mult * i] * c1_5tap + piSrc[mult * i] * c2_5tap +
                     piSrc[mult * i + 1] * c3_5tap + piSrc[mult * i + strOffset] * c4_5tap + offset_5tap) >>
                    shift_5tap;
        } else {
          pDst[i] = (piSrc[mult * i - strOffset] * c0_5tap + piSrc[mult * i] * c1_5tap + piSrc[mult * i - 1] * c2_5tap +
                     piSrc[mult * i + 1] * c3_5tap + piSrc[mult * i + strOffset] * c4_5tap + offset_5tap) >>
                    shift_5tap;
        }
      } else {
        piSrc = pRecSrc0 - iRecStride2;

        if ((i == 0 && !bLeftAvaillable) || (i == avaiAboveSizes - 1 + logSubWidthC)) {
          pDst[i] = ((piSrc[mult * i] * c0_6tap + piSrc[mult * i] * c1_6tap + piSrc[mult * i + 1] * c2_6tap) +
                     (piSrc[mult * i + strOffset] * c3_6tap + piSrc[mult * i + strOffset] * c4_6tap +
                      piSrc[mult * i + 1 + strOffset] * c5_6tap) +
                     offset_6tap) >>
                    shift_6tap;
        } else {
          pDst[i] = ((piSrc[mult * i] * c0_6tap + piSrc[mult * i - 1] * c1_6tap + piSrc[mult * i + 1] * c2_6tap) +
                     (piSrc[mult * i + strOffset] * c3_6tap + piSrc[mult * i - 1 + strOffset] * c4_6tap +
                      piSrc[mult * i + 1 + strOffset] * c5_6tap) +
                     offset_6tap) >>
                    shift_6tap;
        }
      }
    }
  }

  if (bLeftAvaillable) {
    pDst = pDst0 - 1;

    piSrc = pRecSrc0 - 2 - logSubWidthC;

    int availlableLeftSizes = availlableLeftUnit * chromaUnitHeight;
    for (int j = 0; j < availlableLeftSizes; j++) {
      if (pu.cs->sps->getCclmCollocatedChromaFlag()) {
        if ((j == 0 && !bAboveAvaillable) || (j == availlableLeftSizes - 1 + logSubWidthC)) {
          pDst[0] = (piSrc[1] * c0_5tap + piSrc[1] * c1_5tap + piSrc[0] * c2_5tap + piSrc[2] * c3_5tap +
                     piSrc[1 + strOffset] * c4_5tap + offset_5tap) >>
                    shift_5tap;
        } else {
          pDst[0] = (piSrc[1 - strOffset] * c0_5tap + piSrc[1] * c1_5tap + piSrc[0] * c2_5tap + piSrc[2] * c3_5tap +
                     piSrc[1 + strOffset] * c4_5tap + offset_5tap) >>
                    shift_5tap;
        }
      } else {
        pDst[0] = ((piSrc[1] * c0_6tap + piSrc[0] * c1_6tap + piSrc[2] * c2_6tap) +
                   (piSrc[1 + strOffset] * c3_6tap + piSrc[strOffset] * c4_6tap + piSrc[2 + strOffset] * c5_6tap) +
                   offset_6tap) >>
                  shift_6tap;
      }

      piSrc += iRecStride2;
      pDst += iDstStride;
    }
  }

  if (pu.cs->sps->getCclmCollocatedChromaFlag()) {
    // TODO: unroll loop
    for (int j = 0; j < uiCHeight; j++) {
      for (int i = 0; i < uiCWidth; i++) {
        if (i == 0 && !bLeftAvaillable) {
          if (j == 0 && !bAboveAvaillable) {
            pDst0[i] = (pRecSrc0[mult * i] * c0_5tap + pRecSrc0[mult * i] * c1_5tap + pRecSrc0[mult * i] * c2_5tap +
                        pRecSrc0[mult * i + 1] * c3_5tap + pRecSrc0[mult * i + strOffset] * c4_5tap + offset_5tap) >>
                       shift_5tap;
          } else {
            pDst0[i] = (pRecSrc0[mult * i - strOffset] * c0_5tap + pRecSrc0[mult * i] * c1_5tap +
                        pRecSrc0[mult * i] * c2_5tap + pRecSrc0[mult * i + 1] * c3_5tap +
                        pRecSrc0[mult * i + strOffset] * c4_5tap + offset_5tap) >>
                       shift_5tap;
          }
        } else if (j == 0 && !bAboveAvaillable) {
          pDst0[i] = (pRecSrc0[mult * i] * c0_5tap + pRecSrc0[mult * i] * c1_5tap + pRecSrc0[mult * i - 1] * c2_5tap +
                      pRecSrc0[mult * i + 1] * c3_5tap + pRecSrc0[mult * i + strOffset] * c4_5tap + offset_5tap) >>
                     shift_5tap;
        } else {
          pDst0[i] = (pRecSrc0[mult * i - strOffset] * c0_5tap + pRecSrc0[mult * i] * c1_5tap +
                      pRecSrc0[mult * i - 1] * c2_5tap + pRecSrc0[mult * i + 1] * c3_5tap +
                      pRecSrc0[mult * i + strOffset] * c4_5tap + offset_5tap) >>
                     shift_5tap;
        }
      }
      pDst0 += iDstStride;
      pRecSrc0 += iRecStride2;
    }
    return;
  }

#define GET_LUMA_REC_PIX_INC \
  pDst0 += iDstStride;       \
  pRecSrc0 += iRecStride2

#define GET_LUMA_REC_PIX_OP2(ADDR)                                                                                    \
  pDst0[ADDR] =                                                                                                       \
      (pRecSrc0[((ADDR) << logSubWidthC)] * c0_6tap + pRecSrc0[((ADDR) << logSubWidthC) + 1] * c1_6tap +              \
       pRecSrc0[((ADDR) << logSubWidthC) - 1] * c2_6tap + pRecSrc0[((ADDR) << logSubWidthC) + iRecStride] * c3_6tap + \
       pRecSrc0[((ADDR) << logSubWidthC) + 1 + iRecStride] * c4_6tap +                                                \
       pRecSrc0[((ADDR) << logSubWidthC) - 1 + iRecStride] * c5_6tap + offset_6tap) >>                                \
      shift_6tap

#define GET_LUMA_REC_PIX_OP1(ADDR)                                                                                 \
  pDst0[ADDR] =                                                                                                    \
      !(ADDR) ? (pRecSrc0[((ADDR) << logSubWidthC)] * c0_6tap + pRecSrc0[((ADDR) << logSubWidthC) + 1] * c1_6tap + \
                 pRecSrc0[((ADDR) << logSubWidthC)] * c2_6tap +                                                    \
                 pRecSrc0[((ADDR) << logSubWidthC) + iRecStride] * c3_6tap +                                       \
                 pRecSrc0[((ADDR) << logSubWidthC) + 1 + iRecStride] * c4_6tap +                                   \
                 pRecSrc0[((ADDR) << logSubWidthC) + iRecStride] * c5_6tap + offset_6tap) >>                       \
                    shift_6tap                                                                                     \
              : GET_LUMA_REC_PIX_OP2(ADDR)

  int width = uiCWidth;
  int height = uiCHeight;

  GCC_WARNING_DISABLE_sequence_point if (bLeftAvaillable) {
    if (pu.chromaFormat == CHROMA_420) {
      GetLumaRecPixel420(width, height, reinterpret_cast<const Pel *>(pRecSrc0), iRecStride,
                         reinterpret_cast<Pel *>(pDst0), iDstStride);
    } else  // TODO add SIMD for 422,444
    {
      SIZE_AWARE_PER_EL_OP(GET_LUMA_REC_PIX_OP2, GET_LUMA_REC_PIX_INC);
    }
  }
  else {
    SIZE_AWARE_PER_EL_OP(GET_LUMA_REC_PIX_OP1, GET_LUMA_REC_PIX_INC);
  }
  GCC_WARNING_RESET
}

#undef GET_LUMA_REC_PIX_INC
#undef GET_LUMA_REC_PIX_OP1
#undef GET_LUMA_REC_PIX_OP2
#undef SIZE_AWARE_PER_EL_OP

template <typename T>
void IntraPrediction::xGetLMParameters(const PredictionUnit &pu, const ComponentID compID, const CompArea &chromaArea,
                                       int &a, int &b, int &iShift) {
  CHECK(compID == COMPONENT_Y, "");

  const SizeType cWidth = chromaArea.width;
  const SizeType cHeight = chromaArea.height;

  const Position posLT = chromaArea;

  const CodingUnit &cu = pu;
  const CodingStructure &cs = *cu.cs;

  const SPS &sps = *cs.sps;
  const uint32_t tuWidth = chromaArea.width;
  const uint32_t tuHeight = chromaArea.height;
  const ChromaFormat nChromaFormat = sps.getChromaFormatIdc();

  const int baseUnitSize = 1 << MIN_CU_LOG2;
  const int unitWidth = baseUnitSize >> getComponentScaleX(chromaArea.compID, nChromaFormat);
  const int unitHeight = baseUnitSize >> getComponentScaleX(chromaArea.compID, nChromaFormat);

  const int tuWidthInUnits = tuWidth / unitWidth;
  const int tuHeightInUnits = tuHeight / unitHeight;
  const int aboveUnits = tuWidthInUnits;
  const int leftUnits = tuHeightInUnits;
  int topTemplateSampNum = 2 * cWidth;  // for MDLM, the template sample number is 2W or 2H;
  int leftTemplateSampNum = 2 * cHeight;
  CHECKD(!(m_topRefLength >= topTemplateSampNum), "Error!");
  CHECKD(!(m_leftRefLength >= leftTemplateSampNum), "Error!");
  int totalAboveUnits = (topTemplateSampNum + (unitWidth - 1)) / unitWidth;
  int totalLeftUnits = (leftTemplateSampNum + (unitHeight - 1)) / unitHeight;
  int aboveRightUnits = totalAboveUnits - aboveUnits;
  int leftBelowUnits = totalLeftUnits - leftUnits;

  int curChromaMode = pu.intraDir[1];
  bool aboveAvailable = 0, leftAvailable = 0;
  int neighborSize;

  const TransformUnit &tu = *getTU(cu, chromaArea.pos(), CH_C);

  T *srcColor0, *curChroma0;
  int srcStride, curStride;

  T *temp;
  if ((curChromaMode == MDLM_L_IDX) || (curChromaMode == MDLM_T_IDX)) {
    srcStride = 2 * MAX_CU_SIZE + 1;
    temp = (reinterpret_cast<T *>(m_pMdlmTemp)) + srcStride + 1;
  } else {
    srcStride = MAX_CU_SIZE + 1;
    temp = (reinterpret_cast<T *>(m_piTemp)) + srcStride + 1;
  }
  srcColor0 = temp;
  curChroma0 = reinterpret_cast<T *>(getPredictorPtr(compID));

  curStride = m_topRefLength + 1;

  curChroma0 += curStride + 1;

  unsigned internalBitDepth = sps.getBitDepth(CHANNEL_TYPE_CHROMA);

  int minLuma[2] = {MAX_INT, 0};
  int maxLuma[2] = {-MAX_INT, 0};

  T *src = srcColor0 - srcStride;
  T *cur = curChroma0 - curStride;
  int actualTopTemplateSampNum = 0;
  int actualLeftTemplateSampNum = 0;
  if (curChromaMode == MDLM_T_IDX) {
    aboveRightUnits = aboveRightUnits > (cHeight / unitWidth) ? cHeight / unitWidth : aboveRightUnits;
    int avaiAboveUnits =
        isAboveAvailable(tu, CHANNEL_TYPE_CHROMA, posLT, aboveUnits + aboveRightUnits, unitWidth, neighborSize);
    aboveAvailable = avaiAboveUnits >= tuWidthInUnits;
    actualTopTemplateSampNum = unitWidth * avaiAboveUnits;
  } else if (curChromaMode == MDLM_L_IDX) {
    leftBelowUnits = leftBelowUnits > (cWidth / unitHeight) ? cWidth / unitHeight : leftBelowUnits;
    int avaiLeftUnits =
        isLeftAvailable(tu, CHANNEL_TYPE_CHROMA, posLT, leftUnits + leftBelowUnits, unitHeight, neighborSize);
    leftAvailable = avaiLeftUnits >= tuHeightInUnits;
    actualLeftTemplateSampNum = unitHeight * avaiLeftUnits;
  } else if (curChromaMode == LM_CHROMA_IDX) {
    int avaiAboveUnits = isAboveAvailable(tu, CHANNEL_TYPE_CHROMA, posLT, aboveUnits, unitWidth, neighborSize);
    aboveAvailable = avaiAboveUnits == tuWidthInUnits;
    int avaiLeftUnits = isLeftAvailable(tu, CHANNEL_TYPE_CHROMA, posLT, leftUnits, unitHeight, neighborSize);
    leftAvailable = avaiLeftUnits == tuHeightInUnits;
    actualTopTemplateSampNum = cWidth;
    actualLeftTemplateSampNum = cHeight;
  }
  int startPos[2];  // 0:Above, 1: Left
  int pickStep[2];

  int aboveIs4 = leftAvailable ? 0 : 1;
  int leftIs4 = aboveAvailable ? 0 : 1;

  startPos[0] = actualTopTemplateSampNum >> (2 + aboveIs4);
  pickStep[0] = std::max(1, actualTopTemplateSampNum >> (1 + aboveIs4));

  startPos[1] = actualLeftTemplateSampNum >> (2 + leftIs4);
  pickStep[1] = std::max(1, actualLeftTemplateSampNum >> (1 + leftIs4));

  T selectLumaPix[4] = {0, 0, 0, 0};
  T selectChromaPix[4] = {0, 0, 0, 0};

  int cntT, cntL;
  cntT = cntL = 0;
  int cnt = 0;
  if (aboveAvailable) {
    cntT = std::min(actualTopTemplateSampNum, (1 + aboveIs4) << 1);
    src = srcColor0 - srcStride;
    cur = curChroma0 - curStride;
    for (int pos = startPos[0]; cnt < cntT; pos += pickStep[0], cnt++) {
      selectLumaPix[cnt] = src[pos];
      selectChromaPix[cnt] = cur[pos];
    }
  }

  if (leftAvailable) {
    cntL = std::min(actualLeftTemplateSampNum, (1 + leftIs4) << 1);
    src = srcColor0 - 1;
    cur = curChroma0 - 1;
    for (int pos = startPos[1], cnt = 0; cnt < cntL; pos += pickStep[1], cnt++) {
      selectLumaPix[cnt + cntT] = src[pos * srcStride];
      selectChromaPix[cnt + cntT] = cur[pos * curStride];
    }
  }
  cnt = cntL + cntT;

  if (cnt == 2) {
    selectLumaPix[3] = selectLumaPix[0];
    selectChromaPix[3] = selectChromaPix[0];
    selectLumaPix[2] = selectLumaPix[1];
    selectChromaPix[2] = selectChromaPix[1];
    selectLumaPix[0] = selectLumaPix[1];
    selectChromaPix[0] = selectChromaPix[1];
    selectLumaPix[1] = selectLumaPix[3];
    selectChromaPix[1] = selectChromaPix[3];
  }

  int minGrpIdx[2] = {0, 2};
  int maxGrpIdx[2] = {1, 3};
  int *tmpMinGrp = minGrpIdx;
  int *tmpMaxGrp = maxGrpIdx;
  if (selectLumaPix[tmpMinGrp[0]] > selectLumaPix[tmpMinGrp[1]]) std::swap(tmpMinGrp[0], tmpMinGrp[1]);
  if (selectLumaPix[tmpMaxGrp[0]] > selectLumaPix[tmpMaxGrp[1]]) std::swap(tmpMaxGrp[0], tmpMaxGrp[1]);
  if (selectLumaPix[tmpMinGrp[0]] > selectLumaPix[tmpMaxGrp[1]])
    std::swap(tmpMinGrp, tmpMaxGrp);  // TODO: really? not std::swap(tmpMinGrp[0], tmpMaxGrp[1]); ?
  if (selectLumaPix[tmpMinGrp[1]] > selectLumaPix[tmpMaxGrp[0]]) std::swap(tmpMinGrp[1], tmpMaxGrp[0]);

  minLuma[0] = (selectLumaPix[tmpMinGrp[0]] + selectLumaPix[tmpMinGrp[1]] + 1) >> 1;
  minLuma[1] = (selectChromaPix[tmpMinGrp[0]] + selectChromaPix[tmpMinGrp[1]] + 1) >> 1;
  maxLuma[0] = (selectLumaPix[tmpMaxGrp[0]] + selectLumaPix[tmpMaxGrp[1]] + 1) >> 1;
  maxLuma[1] = (selectChromaPix[tmpMaxGrp[0]] + selectChromaPix[tmpMaxGrp[1]] + 1) >> 1;

  if (leftAvailable || aboveAvailable) {
    int diff = maxLuma[0] - minLuma[0];
    if (diff > 0) {
      int diffC = maxLuma[1] - minLuma[1];
      int x = getLog2(diff);
      static const uint8_t DivSigTable[1 << 4] = {// 4bit significands - 8 ( MSB is omitted )
                                                  0, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 0};
      int normDiff = (diff << 4 >> x) & 15;
      int v = DivSigTable[normDiff] | 8;
      x += normDiff != 0;

      int y = getLog2(abs(diffC)) + 1;
      int add = 1 << y >> 1;
      a = (diffC * v + add) >> y;
      iShift = 3 + x - y;
      if (iShift < 1) {
        iShift = 1;
        a = ((a == 0) ? 0 : (a < 0) ? -15 : 15);  // a=Sign(a)*15
      }
      b = minLuma[1] - ((a * minLuma[0]) >> iShift);
    } else {
      a = 0;
      b = minLuma[1];
      iShift = 0;
    }
  } else {
    a = 0;

    b = 1 << (internalBitDepth - 1);

    iShift = 0;
  }
}

void IntraPrediction::initIntraMip(const PredictionUnit &pu, const CompArea &area) {
  CHECK(area.width > MIP_MAX_WIDTH || area.height > MIP_MAX_HEIGHT, "Error: block size not supported for MIP");

  // prepare input (boundary) data for prediction
  Pel *ptrSrc = getPredictorPtr(area.compID);
  const int srcStride = m_topRefLength + 1;  // TODO: check this if correct
  const int srcHStride = m_leftRefLength + 1;

  m_matrixIntraPred.prepareInputForPred(CPelBuf(ptrSrc, srcStride, srcHStride), area,
                                        pu.slice->getSPS()->getBitDepth(toChannelType(area.compID)), area.compID);
}

void IntraPrediction::predIntraMip(const ComponentID compId, PelBuf &piPred, const PredictionUnit &pu) {
  CHECK(piPred.width > MIP_MAX_WIDTH || piPred.height > MIP_MAX_HEIGHT, "Error: block size not supported for MIP");
  CHECK(piPred.width != (1 << getLog2(piPred.width)) || piPred.height != (1 << getLog2(piPred.height)),
        "Error: expecting blocks of size 2^M x 2^N");

  // generate mode-specific prediction
  uint32_t modeIdx = MAX_NUM_MIP_MODE;
  bool transposeFlag = false;
  if (compId == COMPONENT_Y) {
    modeIdx = pu.intraDir[CHANNEL_TYPE_LUMA];
    transposeFlag = pu.mipTransposedFlag();
  } else {
    const PredictionUnit &coLocatedLumaPU = PU::getCoLocatedLumaPU(pu);

    CHECK(pu.intraDir[CHANNEL_TYPE_CHROMA] != DM_CHROMA_IDX, "Error: MIP is only supported for chroma with DM_CHROMA.");
    CHECK(!coLocatedLumaPU.mipFlag(), "Error: Co-located luma CU should use MIP.");

    modeIdx = coLocatedLumaPU.intraDir[CHANNEL_TYPE_LUMA];
    transposeFlag = coLocatedLumaPU.mipTransposedFlag();
  }

  CHECK(modeIdx >= getNumModesMip(piPred), "Error: Wrong MIP mode index");

  const int bitDepth = pu.slice->getSPS()->getBitDepth(toChannelType(compId));
  m_matrixIntraPred.predBlock(piPred, modeIdx, piPred, transposeFlag, bitDepth, compId);
}

//! \}
