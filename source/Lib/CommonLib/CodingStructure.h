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

/** \file     CodingStructure.h
 *  \brief    A class managing the coding information for a specific image part
 */

#ifndef __CODINGSTRUCTURE__
#define __CODINGSTRUCTURE__

#include "Unit.h"
#include "Buffer.h"
#include "CommonDef.h"
#include "UnitPartitioner.h"
#include "Slice.h"

#include <vector>

struct Picture;

enum PictureType {
  PIC_RECONSTRUCTION,
  PIC_RECON_WRAP,
  // PIC_PREDICTION,
  NUM_PIC_TYPES
};

extern ThreadSafeCUCache g_globalUnitCache;

#define NUM_PARTS_IN_CTU (MAX_CU_SIZE * MAX_CU_SIZE) >> (MIN_CU_LOG2 << 1)

struct CtuLoopFilterData {
  LoopFilterParam lfParam[NUM_PARTS_IN_CTU];
};

struct CtuCuPtr {
  CodingUnit* cuPtr[MAX_NUM_CHANNEL_TYPE][NUM_PARTS_IN_CTU];
};

struct CtuData {
  // SAOBlkParam          saoParam;

  // CodingUnit*     cuPtr  [MAX_NUM_CHANNEL_TYPE][NUM_PARTS_IN_CTU];
  // LoopFilterParam lfParam[NUM_EDGE_DIR]        [NUM_PARTS_IN_CTU];
  MotionInfo motion[NUM_PARTS_IN_CTU];
};

struct CtuAlfData {
  uint8_t ccAlfFilterControl[2];
  uint8_t alfCtuEnableFlag[MAX_NUM_COMPONENT];
  uint8_t alfCtuAlternative[MAX_NUM_COMPONENT - 1];
  short alfCtbFilterIndex;
  CtuAlfData() { reset(); }
  void reset() {
    ccAlfFilterControl[0] = ccAlfFilterControl[1] = 0;
    alfCtuEnableFlag[0] = alfCtuEnableFlag[1] = alfCtuEnableFlag[2] = 0;
    alfCtuAlternative[0] = alfCtuAlternative[1] = 0;
    alfCtbFilterIndex = 0;
  }
};

// ---------------------------------------------------------------------------
// coding structure
// ---------------------------------------------------------------------------

class CodingStructure {
 public:
  UnitArea area;

  Picture* picture;

  UnitScale unitScale[MAX_NUM_COMPONENT];
  int chromaQpAdj;
  std::shared_ptr<const VPS> vps;
  std::shared_ptr<const SPS> sps;
  std::shared_ptr<const PPS> pps;
  PicHeader* picHeader;
  std::shared_ptr<APS> alfApss[ALF_CTB_MAX_NUM_APS];
  std::shared_ptr<APS> lmcsAps;
  const PreCalcValues* pcv;
  CtuData* m_ctuData;
  SAOBlkParam* m_saoBlkParam;
  CtuLoopFilterData* m_ctuLoopFilterDataHorEdge;
  CtuCuPtr* m_ctuCuPtr;
  size_t m_ctuDataSize;

  CodingStructure(/*std::shared_ptr<CUCache>, std::shared_ptr<TUCache>*/);

  void create(const UnitArea& _unit);
  void create(const ChromaFormat& _chromaFormat, const Area& _area);
  void destroy();

  void rebindPicBufs();

  // ---------------------------------------------------------------------------
  // global accessors
  // ---------------------------------------------------------------------------

#if _DEBUG
  const CodingUnit* getCU(Position pos, ChannelType _chType) const {
    if (area.blocks[_chType].contains(pos)) {
      ptrdiff_t rsAddr = ctuRsAddr(pos, _chType);
      ptrdiff_t inCtu = inCtuPos(pos, _chType);
      return m_ctuCuPtr[rsAddr].cuPtr[_chType][inCtu];
    } else
      return nullptr;
  }

  CodingUnit* getCU(Position pos, ChannelType _chType) {
    if (area.blocks[_chType].contains(pos)) {
      ptrdiff_t rsAddr = ctuRsAddr(pos, _chType);
      ptrdiff_t inCtu = inCtuPos(pos, _chType);
      return m_ctuCuPtr[rsAddr].cuPtr[_chType][inCtu];
    } else
      return nullptr;
  }
  const CodingUnit* getCUFast(Position pos, ChannelType _chType) const {
    ptrdiff_t rsAddr = ctuRsAddr(pos, _chType);
    ptrdiff_t inCtu = inCtuPos(pos, _chType);
    return m_ctuCuPtr[rsAddr].cuPtr[_chType][inCtu];
  }

  CodingUnit* getCUFast(Position pos, ChannelType _chType) {
    ptrdiff_t rsAddr = ctuRsAddr(pos, _chType);
    ptrdiff_t inCtu = inCtuPos(pos, _chType);
    return m_ctuCuPtr[rsAddr].cuPtr[_chType][inCtu];
  }
#else
  const CodingUnit* getCU(Position pos, ChannelType _chType) const {
    if (area.blocks[_chType].contains(pos))
      return m_ctuCuPtr[ctuRsAddr(pos, _chType)].cuPtr[_chType][inCtuPos(pos, _chType)];
    else
      return nullptr;
  }
  CodingUnit* getCU(Position pos, ChannelType _chType) {
    if (area.blocks[_chType].contains(pos))
      return m_ctuCuPtr[ctuRsAddr(pos, _chType)].cuPtr[_chType][inCtuPos(pos, _chType)];
    else
      return nullptr;
  }

  const CodingUnit* getCUFast(Position pos, ChannelType _chType) const {
    return m_ctuCuPtr[ctuRsAddr(pos, _chType)].cuPtr[_chType][inCtuPos(pos, _chType)];
  }
  CodingUnit* getCUFast(Position pos, ChannelType _chType) {
    return m_ctuCuPtr[ctuRsAddr(pos, _chType)].cuPtr[_chType][inCtuPos(pos, _chType)];
  }
#endif

  const CodingUnit* getCURestricted(const Position& pos, const Position curPos, const unsigned curSliceIdx,
                                    const unsigned curTileIdx, const ChannelType _chType) const;
  const CodingUnit* getCURestricted(const Position& pos, const CodingUnit& curCu, const ChannelType _chType,
                                    const CodingUnit* guess = nullptr) const;

  CodingUnit& addCU(const UnitArea& unit, const ChannelType _chType, const TreeType treeType, const ModeType modeType,
                    const CodingUnit* cuLeft, const CodingUnit* cuAbove);
  TransformUnit& addTU(const UnitArea& unit, const ChannelType _chType, CodingUnit& cu);
  void addEmptyTUs(Partitioner& partitioner, CodingUnit& cu);
  CUTraverser traverseCUs(const UnitArea& _unit);

  cCUTraverser traverseCUs(const UnitArea& _unit) const;

  void initStructData();

 private:
  void createInternals(const UnitArea& _unit);

 private:
  // Pel*            m_predBuf    [MAX_NUM_COMPONENT];

  Mv* m_dmvrMvCache;
  size_t m_dmvrMvCacheSize;
  ptrdiff_t m_dmvrMvCacheOffset;

  unsigned m_numCUs;
  unsigned m_numTUs;

  std::shared_ptr<CUCache> m_cuCache;
  std::shared_ptr<TUCache> m_tuCache;

  PelStorage m_reco;
  PelStorage m_rec_wrap;
  // PelStorage m_pred;

  CodingUnit* m_lastCU = nullptr;

  size_t m_widthInCtus;
  PosType m_ctuSizeMask[2];
  PosType m_ctuWidthLog2[2];
  Pel* m_pictureResidualBuffer;
  ptrdiff_t m_residualOffset;

 public:
  void setupCUTUCache(std::shared_ptr<CUCache> cuCache, std::shared_ptr<TUCache> tuCache) {
    m_cuCache = cuCache;
    m_tuCache = tuCache;
  }

  void releaseCUTUCache() {
    m_cuCache.reset();
    m_tuCache.reset();
  }

  // in CTU coordinates
  int ctuRsAddr(int col, int line) const { return col + (line * (int)m_widthInCtus); }
  // in sample coordinates
  int ctuRsAddr(Position pos, ChannelType chType) const {
    Position posL = recalcPosition(area.chromaFormat, chType, CH_L, pos);
    return ctuRsAddr(posL.x >> pcv->maxCUWidthLog2, posL.y >> pcv->maxCUHeightLog2);
  }
  // 4x4 luma block position within the CTU
  int inCtuPos(Position pos, ChannelType chType) const {
    return (unitScale[chType].scaleHor(pos.x) & m_ctuSizeMask[chType]) +
           ((unitScale[chType].scaleVer(pos.y) & m_ctuSizeMask[chType]) << m_ctuWidthLog2[chType]);
  }

  CtuData& getCtuData(int col, int line) { return m_ctuData[ctuRsAddr(col, line)]; }
  const CtuData& getCtuData(int col, int line) const { return m_ctuData[ctuRsAddr(col, line)]; }
  CtuLoopFilterData& getCtuLoopFilterDataHorEdge(int col, int line) {
    return m_ctuLoopFilterDataHorEdge[ctuRsAddr(col, line)];
  }
  const CtuLoopFilterData& getCtuLoopFilterDataHorEdge(int col, int line) const {
    return m_ctuLoopFilterDataHorEdge[ctuRsAddr(col, line)];
  }

  CtuData& getCtuData(int addr) { return m_ctuData[addr]; }
  const CtuData& getCtuData(int addr) const { return m_ctuData[addr]; }
  CtuCuPtr& getCtuCuPtrData(int addr) { return m_ctuCuPtr[addr]; }
  const CtuCuPtr& getCtuCuPtrData(int addr) const { return m_ctuCuPtr[addr]; }
  CtuCuPtr& getCtuCuPtrData(int col, int line) { return m_ctuCuPtr[ctuRsAddr(col, line)]; }
  const CtuCuPtr& getCtuCuPtrData(int col, int line) const { return m_ctuCuPtr[ctuRsAddr(col, line)]; }

  int m_IBCBufferWidth;
  std::vector<PelStorage> m_virtualIBCbuf;
  void initVIbcBuf(int numCtuLines, ChromaFormat chromaFormatIDC, int ctuSize);
  void fillIBCbuffer(CodingUnit& cu, int lineIdx);

  // PelStorage *m_alfBuf;

  MotionBuf getMotionBuf(const Area& _area);
  MotionBuf getMotionBuf(const UnitArea& _area) { return getMotionBuf(_area.Y()); }

  const CMotionBuf getMotionBuf(const Area& _area) const;
  const CMotionBuf getMotionBuf(const UnitArea& _area) const { return getMotionBuf(_area.Y()); }

  MotionInfo& getMotionInfo(const Position& pos) {
    return getCtuData(ctuRsAddr(pos, CH_L)).motion[inCtuPos(pos, CH_L)];
  }
  const MotionInfo& getMotionInfo(const Position& pos) const {
    return getCtuData(ctuRsAddr(pos, CH_L)).motion[inCtuPos(pos, CH_L)];
  }

  LoopFilterParam const* getLFPMapPtrHorEdge(ptrdiff_t _ctuRsAddr) const {
    return m_ctuLoopFilterDataHorEdge[_ctuRsAddr].lfParam;
  }
  LoopFilterParam* getLFPMapPtrHorEdge(ptrdiff_t _ctuRsAddr) { return m_ctuLoopFilterDataHorEdge[_ctuRsAddr].lfParam; }
  LoopFilterParam* getLFPMapPtrHorEdge(const Position& pos, ChannelType chType) {
    return getLFPMapPtrHorEdge(ctuRsAddr(pos, chType)) + inCtuPos(pos, chType);
  }
  ptrdiff_t getLFPMapStride() const { return (ptrdiff_t(1) << m_ctuWidthLog2[CH_L]); }

  UnitScale getScaling(const UnitScale::ScaliningType type, const ChannelType chType = CH_L) const {
    return type == UnitScale::MI_MAP ? g_miScaling : unitScale[chType];
  }

 public:
#if ADAPTIVE_BIT_DEPTH
  PelBuf getRecoBuf(const CompArea& blk, int bytePerPixel) {
    return m_reco.bufs[blk.compID].subBuf(blk.pos(), blk.size(), bytePerPixel);
  }

  PelUnitBuf getRecoBuf(const UnitArea& unit, int bytePerPixel) { return m_reco.subBuf(unit, bytePerPixel); }
#endif
  PelBuf getRecoBuf(const CompArea& blk) { return m_reco.bufs[blk.compID].subBuf(blk.pos(), blk.size()); }
  const CPelBuf getRecoBuf(const CompArea& blk) const { return m_reco.bufs[blk.compID].subBuf(blk.pos(), blk.size()); }
  PelUnitBuf getRecoBuf(const UnitArea& unit) { return m_reco.subBuf(unit); }
  const CPelUnitBuf getRecoBuf(const UnitArea& unit) const { return m_reco.subBuf(unit); }

  void allocateTUResidualBuffer(TransformUnit& tu, int compID) {
    size_t area = tu.blocks[compID].area();
    tu.residuals[compID] = m_pictureResidualBuffer + m_residualOffset;
    m_residualOffset += area;
  }
  void setupPictureResidualBuffer(Pel* buffer) {
    m_pictureResidualBuffer = buffer;
    m_residualOffset = 0;
  }
  Pel* getPictureResidualBuffer() { return m_pictureResidualBuffer; }
  size_t getResidualPictureSize() {
    int w = pcv->lumaWidth;
    int h = pcv->lumaHeight;
    size_t s = w * h;
    switch (pcv->chrFormat) {
      case CHROMA_420:
        s += s / 2;
        break;
      case CHROMA_422:
        s += s;
        break;
      case CHROMA_444:
        s = s * 3;
      default:
        break;
    }
    return s + 16;
  }

  PelStorage* m_alfStorage;
  PelStorage* getAlfDstPtr() { return m_alfStorage; }
  void setAlfDstPtr(PelStorage* ptr) { m_alfStorage = ptr; }

  // reco buffer
  PelBuf getRecoBuf(const ComponentID compID, bool wrap = false) {
    return wrap ? m_rec_wrap.get(compID) : m_reco.get(compID);
  }
  const CPelBuf getRecoBuf(const ComponentID compID, bool wrap = false) const {
    return wrap ? m_rec_wrap.get(compID) : m_reco.get(compID);
  }
  PelUnitBuf getRecoBuf(bool wrap = false) { return wrap ? m_rec_wrap : m_reco; }
  const CPelUnitBuf getRecoBuf(bool wrap = false) const { return wrap ? m_rec_wrap : m_reco; }

  PelStorage* getReconBufPtr() { return &m_reco; }

  PelUnitBuf getPredBuf(const CodingUnit& cu);
  const CPelUnitBuf getPredBuf(const CodingUnit& cu) const;
};

#endif
