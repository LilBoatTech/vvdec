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

/** \file     LoopFilter.h
    \brief    deblocking filter (header)
*/

#ifndef __LOOPFILTER__
#define __LOOPFILTER__

#include "CommonDef.h"
#include "Unit.h"
#include "Picture.h"

//! \ingroup CommonLib
//! \{

#define DEBLOCK_SMALLEST_BLOCK 8

// ====================================================================================================================
// Class definition
// ====================================================================================================================
typedef int (*pfunc_xCalBsSameRef)(const MotionInfo& miP, const MotionInfo& miQ, const Picture* piRefP0,
                                   const Picture* piRefP1, const Picture* piRefQ0, const Picture* piRefQ1);
typedef void (*pfunc_xPelFilterLuma)(Pel* piSrc, const ptrdiff_t step, const ptrdiff_t offset, const int tc,
                                     const bool sw, const int iThrCut, const bool bFilterSecondP,
                                     const bool bFilterSecondQ, const ClpRng& clpRng);
typedef void (*pfunc_xFilteringPandQ)(Pel* src, ptrdiff_t step, const ptrdiff_t offset, int numberPSide,
                                      int numberQSide, int tc);
/// deblocking filter class
class LoopFilter {
 private:
  /// CU-level deblocking function
  template <DeblockEdgeDir edgeDir>
  void xDeblockCtuArea(CodingStructure& cs, const UnitArea& area, const ChannelType chType,
                       LoopFilterParam* loopFilterParamVerEdge);

  // set / get functions
  LFCUParam xGetLoopfilterParam(const CodingUnit& cu);

  // filtering functions
  template <DeblockEdgeDir edgeDir>
  void xGetBoundaryStrengthSingle(LoopFilterParam& lfp, const CodingUnit& cu, const Position& localPos,
                                  const CodingUnit& cuP, CtuData& ctuData, bool pqSameCtu);
  template <DeblockEdgeDir edgeDir>
  void xSetEdgeFilterInsidePu(const CodingUnit& cu, const Area& area, const bool bValue, CtuData& ctuData,
                              LoopFilterParam* loopFilterParamVerEdge);

  template <DeblockEdgeDir edgeDir>
  void xSetMaxFilterLengthPQFromTransformSizes(const CodingUnit& cu, const TransformUnit& currTU, const bool bValue,
                                               bool deriveBdStrngt, CtuData& ctuData, bool pqSameCtu,
                                               LoopFilterParam* loopFilterParamVerEdge);
  template <DeblockEdgeDir edgeDir>
  void xSetMaxFilterLengthPQForCodingSubBlocks(const CodingUnit& cu, CtuData& ctuData,
                                               LoopFilterParam* loopFilterParamVerEdge);

  template <DeblockEdgeDir edgeDir, typename T>
  void xEdgeFilterLuma(CodingStructure& cs, const Position& pos, const LoopFilterParam& lfp);
  template <DeblockEdgeDir edgeDir, typename T>
  void xEdgeFilterChroma(CodingStructure& cs, const Position& pos, const LoopFilterParam& lfp);

#if LUMA_ADAPTIVE_DEBLOCKING_FILTER_QP_OFFSET
  template <typename T>
  void deriveLADFShift(const T* src, const ptrdiff_t stride, int& shift, const DeblockEdgeDir edgeDir, const SPS sps);

#endif
  static const uint16_t sm_tcTable[MAX_QP + 3];
  static const uint8_t sm_betaTable[MAX_QP + 1];

  int (*xCalBsSameRef)(const MotionInfo& miP, const MotionInfo& miQ, const Picture* piRefP0, const Picture* piRefP1,
                       const Picture* piRefQ0, const Picture* piRefQ1);

  void (*xEdgeFilterLumaImp)(Pel* _piSrc, ptrdiff_t srcStep, ptrdiff_t offset, bool sidePisLarge, bool sideQisLarge,
                             int iBeta, int iTc, int maxFilterLengthP, int maxFilterLengthQ, int iSideThreshold,
                             int iThrCut, const ClpRng& clpRng);

  inline bool isCrossedByVirtualBoundariesHorEdge(const PicHeader* picHeader, int y, int height, int& numHorVirBndry,
                                                  int horVirBndryPos[]);
  inline bool isCrossedByVirtualBoundariesVerEdge(const PicHeader* picHeader, int x, int width, int& numVerVirBndry,
                                                  int verVirBndryPos[]);

  inline void xDeriveEdgefilterParam(const int pos, const int numVirBndry, const int virBndryPos[], bool& edgeFilter);

 public:
  LoopFilter();
  ~LoopFilter();
  /// picture-level deblocking filter
  void loopFilterCTU(CodingStructure& cs, const ChannelType chType, const int ctuCol, const int ctuLine,
                     DeblockEdgeDir edgeDir = NUM_EDGE_DIR, LoopFilterParam* loopFilterParamVerEdge = nullptr);

  void calcFilterStrengthsHorEdge(const CodingUnit& cu);
  void calcFilterStrengthsVerEdge(const CodingUnit& cu, LoopFilterParam* loopFilterParamVerEdge);
  void calcFilterStrengthsVerEdgeCTU(CodingStructure& cs, const UnitArea& ctuArea,
                                     LoopFilterParam* loopFilterParamVerEdge);

  int getBeta(const int qp) {
    const int indexB = Clip3(0, MAX_QP, qp);
    return sm_betaTable[indexB];
  }

  void initLoopFilter(int bitDepth);
#ifdef TARGET_SIMD_X86
  void initLoopFilterX86(int bitDepth);
  template <X86_VEXT vext>
  void _initLoopFilterX86(int bitDepth);
#endif
};

//! \}

#endif
