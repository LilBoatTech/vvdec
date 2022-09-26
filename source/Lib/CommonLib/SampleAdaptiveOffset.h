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

/** \file     SampleAdaptiveOffset.h
    \brief    sample adaptive offset class (header)
*/

#ifndef __SAMPLEADAPTIVEOFFSET__
#define __SAMPLEADAPTIVEOFFSET__

#include "CommonDef.h"
#include "Unit.h"
#include "Reshape.h"

#include <array>

//! \ingroup CommonLib
//! \{

template <typename T>
static inline int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// ====================================================================================================================
// Constants
// ====================================================================================================================

#define MAX_SAO_TRUNCATED_BITDEPTH 10

// ====================================================================================================================
// Class definition
// ====================================================================================================================

class SampleAdaptiveOffset {
 public:
  SampleAdaptiveOffset() {
    m_numberOfComponents = 3;  // at most 3
    for (int compIdx = 0; compIdx < m_numberOfComponents; compIdx++) {
      m_lineTopBuf[compIdx] = nullptr;
      m_lineLeftBuf[compIdx] = nullptr;
    }
    m_picH = 0;
    m_picW = 0;
    m_chromaFormat = 0;
  }
  ~SampleAdaptiveOffset();

  void create(int picWidth, int picHeight, ChromaFormat format, uint32_t maxCUWidth, uint32_t maxCUHeight,
              uint32_t maxCUDepth);
  void destroy();

  //  void SAOProcess( CodingStructure& cs );
  void SAOProcessCTU(CodingStructure& cs, const UnitArea& ctuArea);

  void SAOPrepareCTULine(CodingStructure& cs, const UnitArea& lineArea);

  static int getMaxOffsetQVal(const int channelBitDepth) {
    return (1 << (std::min<int>(channelBitDepth, MAX_SAO_TRUNCATED_BITDEPTH) - 5)) - 1;
  }  // Table 9-32, inclusive
  void setReshaper(Reshape* p) { m_pcReshape = p; }

  void initSampleAdaptiveOffset(int bitDepth);

#ifdef TARGET_SIMD_X86
  void initSampleAdaptiveOffsetX86(int bitDepth);
  template <X86_VEXT vext>
  void _initSampleAdaptiveOffsetX86(int bitDepth);
#endif

 protected:
  void deriveLoopFilterBoundaryAvailibility(CodingStructure& cs, const Position& pos, bool& isLeftAvail,
                                            bool& isRightAvail, bool& isAboveAvail, bool& isBelowAvail,
                                            bool& isAboveLeftAvail, bool& isAboveRightAvail, bool& isBelowLeftAvail,
                                            bool& isBelowRightAvail) const;

  template <typename T>
  static void offsetBlock_core(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset, int startIdx,
                               T* resBlk, ptrdiff_t resStride, int width, int height, bool isLeftAvail,
                               bool isRightAvail, bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail,
                               bool isAboveRightAvail, bool isBelowLeftAvail, bool isBelowRightAvail,
                               bool isCtuCrossedByVirtualBoundaries, int horVirBndryPos[], int verVirBndryPos[],
                               int numHorVirBndry, int numVerVirBndry, const T* topline, const T* leftline);

  void offsetBlock_core(const int channelBitDepth, const ClpRng& clpRng, int typeIdx, int* offset, int startIdx,
                        Pel* resBlk, ptrdiff_t resStride, int width, int height, bool isLeftAvail, bool isRightAvail,
                        bool isAboveAvail, bool isBelowAvail, bool isAboveLeftAvail, bool isAboveRightAvail,
                        bool isBelowLeftAvail, bool isBelowRightAvail, const Pel* topline, const Pel* leftline);
  void (*processBO)(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset, int startIdx,
                    const ClpRng& clpRng, int channelBitDepth);
  void (*processEO0)(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset, const ClpRng& clpRng,
                     bool isLeftAvail, bool isRightAvail, const Pel* leftLine);
  void (*processEO90)(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset, const ClpRng& clpRng,
                      bool isAboveAvail, bool isBottomAvail, const Pel* topLine);
  void (*processEO135)(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset,
                       const ClpRng& clpRng, bool isLeftAvail, bool isRightAvail, bool isAboveLeftAvail,
                       bool isAboveAvail, bool isBelowAvail, bool isBelowRightAvail, const Pel* topLine,
                       const Pel* leftLine);
  void (*processEO45)(Pel* resLine, ptrdiff_t resStride, int width, int height, const int* offset, const ClpRng& clpRng,
                      bool isLeftAvail, bool isRightAvail, bool isAboveRightAvail, bool isAboveAvail, bool isBelowAvail,
                      bool isBelowLeftAvail, const Pel* topLine, const Pel* leftLine);

  // void invertQuantOffsets     ( ComponentID compIdx, int typeIdc, int typeAuxInfo, int* dstOffsets, int* srcOffsets )
  // const;
  void reconstructBlkSAOParam(SAOBlkParam& recParam, SAOBlkParam* mergeList[]) const;
  static int getMergeList(CodingStructure& cs, int ctuRsAddr, SAOBlkParam* mergeList[NUM_SAO_MERGE_TYPES]);
  void offsetCTU(const UnitArea& area, PelUnitBuf& res, int ctuRSAddr, CodingStructure& cs);

  template <typename T>
  void saveCTULine(CodingStructure& cs, const UnitArea& area, int compIdx, bool bSaveBottom, bool bSaveRight,
                   bool bSaveRightFast);

  static bool isCrossedByVirtualBoundaries(const PicHeader* picHeader, const Area& area, int& numHorVirBndry,
                                           int& numVerVirBndry, int horVirBndryPos[], int verVirBndryPos[]);
  static bool isProcessDisabled(int xPos, int yPos, int numVerVirBndry, int numHorVirBndry, int verVirBndryPos[],
                                int horVirBndryPos[]);
  Reshape* m_pcReshape;

 protected:
  // std::array<uint32_t, MAX_NUM_COMPONENT> m_offsetStepLog2;   // offset step
  // PelStorage                              m_tempBuf;
  uint32_t m_numberOfComponents;
  uint32_t m_picW;
  uint32_t m_picH;
  uint32_t m_chromaFormat;
  // std::vector<int8_t> signLineBuf1;
  // std::vector<int8_t> signLineBuf2;
  Pel** m_lineLeftBuf[3];
  Pel** m_lineTopBuf[3];
};

//! \}
#endif
