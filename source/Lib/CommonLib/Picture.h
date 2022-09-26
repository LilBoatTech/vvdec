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

/** \file     Picture.h
 *  \brief    Description of a coded picture
 */

#ifndef __PICTURE__
#define __PICTURE__
#include <thread>
#include "Utilities/threading.h"
#include "CommonDef.h"

#include "Common.h"
#include "Unit.h"
#include "Buffer.h"
#include "Unit.h"
#include "Slice.h"
#include "CodingStructure.h"
#include <deque>

#include "InterpolationFilter.h"

class SEI;
typedef std::list<SEI*> SEIMessages;

struct Picture : public UnitArea {
  Picture() = default;
  ~Picture() = default;
#if ADAPTIVE_BIT_DEPTH
  void create(const ChromaFormat& _chromaFormat, const Size& size, const unsigned _maxCUSize, const unsigned margin,
              const int layerId, bool bEnableWrapAround, int bytePerPixel);
#else
  void create(const ChromaFormat& _chromaFormat, const Size& size, const unsigned _maxCUSize, const unsigned margin,
              const int layerId, bool bEnableWrapAround);
#endif
  void resetForUse();
  void destroy();

  // get reference pel pointer from
  const Pel* getRecoBufPtr(const ComponentID compID, bool wrap = false) const { return m_finalBuffer.bufs[compID].buf; }
  const ptrdiff_t getRecoBufStride(const ComponentID compID, bool wrap = false) const {
    return m_finalBuffer.bufs[compID].stride;
  }
  CPelBuf getReconBufFinal(const ComponentID compID, bool wrap = false) const { return m_finalBuffer.bufs[compID]; }

#if ADAPTIVE_BIT_DEPTH
  PelBuf getRecoBuf(const CompArea& blk, bool wrap, int bytePerPixel);
#endif

  PelBuf getRecoBuf(const ComponentID compID, bool wrap = false);
  const CPelBuf getRecoBuf(const ComponentID compID, bool wrap = false) const;
  PelBuf getRecoBuf(const CompArea& blk, bool wrap = false);
  const CPelBuf getRecoBuf(const CompArea& blk, bool wrap = false) const;
  PelUnitBuf getRecoBuf(const UnitArea& unit, bool wrap = false);
  const CPelUnitBuf getRecoBuf(const UnitArea& unit, bool wrap = false) const;
  PelUnitBuf getRecoBuf(bool wrap = false);
  const CPelUnitBuf getRecoBuf(bool wrap = false) const;

  PelBuf getBuf(const ComponentID compID, const PictureType& type) { return m_bufs[type].bufs[compID]; }
  const CPelBuf getBuf(const ComponentID compID, const PictureType& type) const { return m_bufs[type].bufs[compID]; }
  PelBuf getBuf(const CompArea& blk, const PictureType& type);
  const CPelBuf getBuf(const CompArea& blk, const PictureType& type) const;
  PelUnitBuf getBuf(const UnitArea& unit, const PictureType& type);
  const CPelUnitBuf getBuf(const UnitArea& unit, const PictureType& type) const;

  void extendRowBorder(int row);
  void extendPicBorder(bool top = true, bool bottom = true, bool leftrightT = true, bool leftrightB = true,
                       ChannelType chType = MAX_NUM_CHANNEL_TYPE);
  template <typename T>
  void extendPicBorderWrap(ComponentID compID, bool top, bool bottom, bool leftrightT, bool leftrightB);
  void (*paddPicBorderBot)(Pel* pi, ptrdiff_t stride, int width, int xmargin, int ymargin);
  void (*paddPicBorderTop)(Pel* pi, ptrdiff_t stride, int width, int xmargin, int ymargin);
  void (*paddPicBorderLeftRight)(Pel* pi, ptrdiff_t stride, int width, int xmargin, int height);
  template <typename T>
  void paddPicBorderLeftRightWrap(Pel* pi, ptrdiff_t stride, int width, int xmargin, int height, int offset);

  void finalInit(const SPS* sps, const PPS* pps, PicHeader* picHeader, APS* alfApss[ALF_CTB_MAX_NUM_APS], APS* lmcsAps,
                 APS* scalingListAps);
  int getPOC() const { return poc; }
  uint64_t getCts() const { return cts; }
  uint64_t getDts() const { return dts; }
  uint32_t getTLayer() const { return layer; }
  uint64_t getNaluBits() const { return bits; }
  bool getRap() const { return rap; }

  void setBorderExtension(bool bFlag) { isBorderExtended = bFlag; }
  Pel* getOrigin(const PictureType& type, const ComponentID compID) const;
  PelBuf getOriginBuf(const PictureType& type, const ComponentID compID);

  int getDecodingOrderNumber() const { return decodingOrderNumber; }
  void setDecodingOrderNumber(const int val) { decodingOrderNumber = val; }

  int getSpliceIdx(uint32_t idx) const { return m_spliceIdx[idx]; }
  void setSpliceIdx(uint32_t idx, int poc) { m_spliceIdx[idx] = poc; }
  void createSpliceIdx(int nums);
  bool getSpliceFull();
  static void sampleRateConv(const Pel* orgSrc, SizeType orgWidth, SizeType orgHeight, ptrdiff_t orgStride,
                             Pel* scaledSrc, SizeType scaledWidth, SizeType scaledHeight, SizeType paddedWidth,
                             SizeType paddedHeight, ptrdiff_t scaledStride, const int bitDepth,
                             const bool useLumaFilter, const bool downsampling = false);

  static void rescalePicture(const CPelUnitBuf& beforeScaling, const Window& confBefore, const PelUnitBuf& afterScaling,
                             const Window& confAfter, const ChromaFormat chromaFormatIDC, const BitDepths& bitDepths,
                             const bool useLumaFilter, const bool downsampling = false);

 public:
  bool m_isSubPicBorderSaved = false;

  PelStorage m_bufSubPicAbove;
  PelStorage m_bufSubPicBelow;
  PelStorage m_bufSubPicLeft;
  PelStorage m_bufSubPicRight;

  void saveSubPicBorder(int POC, int subPicX0, int subPicY0, int subPicWidth, int subPicHeight);
  void extendSubPicBorder(int POC, int subPicX0, int subPicY0, int subPicWidth, int subPicHeight);
  void restoreSubPicBorder(int POC, int subPicX0, int subPicY0, int subPicWidth, int subPicHeight);

  bool getSubPicSaved() { return m_isSubPicBorderSaved; }
  void setSubPicSaved(bool bVal) { m_isSubPicBorderSaved = bVal; }

  void startProcessingTimer();
  void stopProcessingTimer();
  void resetProcessingTime() { m_dProcessingTime = 0; }
  double getProcessingTime() const { return m_dProcessingTime; }

  std::chrono::time_point<std::chrono::steady_clock> m_processingStartTime;
  double m_dProcessingTime = 0;

  bool isBorderExtended = false;
  bool wrapAroundValid = false;
  unsigned wrapAroundOffset = 0;
  bool referenced = false;
  bool reconstructed = false;
  bool inProgress = false;
  bool neededForOutput = false;
  bool wasLost = false;
  bool longTerm = false;
  bool topField = false;
  bool fieldPic = false;
  int skippedDecCount = 0;
  int m_prevQP[MAX_NUM_CHANNEL_TYPE] = {-1, -1};
  bool nonReferencePictureFlag = false;

  // As long as this field is true, the picture will not be reused or deleted.
  // An external application needs to call DecLib::releasePicture(), when it is done using the picture buffer.
  bool lockedByApplication = false;

  int poc = 0;
  uint64_t cts = 0;  // composition time stamp
  uint64_t dts = 0;  // decoding time stamp
  uint32_t layer = std::numeric_limits<uint32_t>::max();
  uint32_t depth = 0;
  int layerId = NOT_VALID;
  NalUnitType eNalUnitType = NAL_UNIT_INVALID;
  uint64_t bits = 0;  // input nal bit count
  bool rap = 0;       // random access point flag
  int decodingOrderNumber = 0;
  std::vector<int> sliceSubpicIdx;
  std::vector<SubPic> subPictures;
  int numSlices = 1;

  bool subLayerNonReferencePictureDueToSTSA = 0;

  int* m_spliceIdx = nullptr;
  int m_ctuNums = 0;

  PelStorage m_bufs[NUM_PIC_TYPES];
  uint32_t margin = 0;
  const Picture* unscaledPic;

  // WaitCounter m_ctuTaskCounter;
  // WaitCounter m_dmvrTaskCounter;
  // WaitCounter m_borderExtTaskCounter;
  // BlockingBarrier done;
#if RECO_WHILE_PARSE
  // Barrier* ctuParsedBarrier = nullptr;
#endif
#if ALLOW_MIDER_LF_DURING_PICEXT
  // CBarrierVec refPicExtDepBarriers;
#endif

    vvdec::VVDecThreadCounter* m_rowReconCounter = nullptr;
    vvdec::VVDecThreadCounter m_rowCompleteCount;
    vvdec::VVDecThreadCounter* m_rowMotionInfoCounter = nullptr;
  PelStorage m_finalBuffer;
  // Barrier parseDone;

  CodingStructure* cs = nullptr;
  std::vector<Slice*> slices;
  SEIMessages SEIs;

  bool isRefScaled(const PPS* pps) const {
    const PPS* pps0 = slices[0]->getPPS();
    const Size& recoSize = m_bufs[PIC_RECONSTRUCTION].bufs[COMPONENT_Y];
    return (recoSize.width != pps->getPicWidthInLumaSamples() || recoSize.height != pps->getPicHeightInLumaSamples()) ||
           ((pps0->getScalingWindow().getWindowEnabledFlag() || pps->getScalingWindow().getWindowEnabledFlag()) &&
            (pps0->getScalingWindow().getWindowLeftOffset() != pps->getScalingWindow().getWindowLeftOffset() ||
             pps0->getScalingWindow().getWindowRightOffset() != pps->getScalingWindow().getWindowRightOffset() ||
             pps0->getScalingWindow().getWindowTopOffset() != pps->getScalingWindow().getWindowTopOffset() ||
             pps0->getScalingWindow().getWindowBottomOffset() != pps->getScalingWindow().getWindowBottomOffset()));
  }

  std::shared_ptr<PicHeader> picHeader;
  void setPicHead(const std::shared_ptr<PicHeader>& ph);

  bool isWrapAroundEnabled(const PPS* pps) const { return pps->getUseWrapAround() && !isRefScaled(pps); }

  void allocateNewSlice();
  Slice* swapSliceObject(Slice* p, uint32_t i);
  void clearSliceBuffer();

#if TRACE_ENABLE_ITT
  __itt_domain* m_itt_decLibInst;
#endif

 public:
  std::vector<CtuAlfData> m_ctuAlfData;
  CtuAlfData& getCtuAlfData(int ctuAddr) { return m_ctuAlfData[ctuAddr]; }
  const CtuAlfData& getCtuAlfData(int ctuAddr) const { return m_ctuAlfData[ctuAddr]; }
  void resizeCtuAlfData(int numEntries) { m_ctuAlfData.resize(numEntries); }

  void initPicture(int bitDepth);
#if ENABLE_SIMD_OPT_PICTURE
  void initPictureX86(int bitDepth);
  template <X86_VEXT vext>
  void _initPictureX86(int bitDepth);
#endif

#if ADAPTIVE_BIT_DEPTH
  int m_bytePerPixel;
#endif
};

int calcAndPrintHashStatus(const CPelUnitBuf& pic, const class SEIDecodedPictureHash* pictureHashSEI,
                           const BitDepths& bitDepths, const MsgLevel msgl);

#endif
