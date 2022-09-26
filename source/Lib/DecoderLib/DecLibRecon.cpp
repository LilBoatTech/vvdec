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

/** \file     DecLibRecon.cpp
    \brief    decoder class
*/
#include <atomic>
#include "DecLib.h"
#include "Utilities/threadpool.h"
#include "CommonLib/TrQuant.h"
#include "CommonLib/TrQuant_EMT.h"
#include "CommonLib/InterPrediction.h"
#include "CommonLib/IntraPrediction.h"
#include "CommonLib/Unit.h"
#include "CommonLib/Buffer.h"
#include "CommonLib/UnitTools.h"

#include "CommonLib/dtrace_next.h"
#include "CommonLib/dtrace_buffer.h"

#ifdef TRACE_ENABLE_ITT
extern __itt_domain* itt_domain_dec;
extern std::vector<__itt_domain*> itt_domain_decInst;

extern __itt_string_handle* itt_handle_alf;
extern __itt_string_handle* itt_handle_presao;
extern __itt_string_handle* itt_handle_sao;
extern __itt_string_handle* itt_handle_lfl;
extern __itt_string_handle* itt_handle_intra;
extern __itt_string_handle* itt_handle_inter;
extern __itt_string_handle* itt_handle_mider;
extern __itt_string_handle* itt_handle_lfcl;
extern __itt_string_handle* itt_handle_ext;
extern __itt_string_handle* itt_handle_dmvr;
extern __itt_string_handle* itt_handle_rsp;

extern __itt_string_handle* itt_handle_schedTasks;
extern __itt_string_handle* itt_handle_waitTasks;

// create global domain for DecLib
extern __itt_domain* itt_domain_glb;
// create a global counter
extern __itt_counter itt_frame_counter;

#  define ITT_TASKSTART(d, t) __itt_task_begin((d), __itt_null, __itt_null, (t))
#  define ITT_TASKEND(d, t) __itt_task_end((d))
#else
#  define ITT_TASKSTART(d, t)
#  define ITT_TASKEND(d, t)
#endif

//! \ingroup DecoderLib
//! \{
DecLibRecon::DecLibRecon() {
    m_rowDependencyBitmap = nullptr;
    m_frameParseRowBitmap = nullptr;
    m_numWords = 0;
    m_numRows = 0;
  m_currFuncBitDepth = -1;

#if ENABLE_SIMD_OPT_BUFFER
  g_pelBufOP.initPelBufOpsX86();
#endif
}

void DecLibRecon::decodeFrame() {
  m_sPerLineMiHist.resize(m_currDecompPic->cs->pcv->heightInCtus);

  // here to performa to decode a frame
  int lines = m_currDecompPic->cs->pcv->heightInCtus;
  int w = m_currDecompPic->cs->pcv->widthInCtus;
  CodingStructure& cs = *m_currDecompPic->cs;
  if (cs.sps->getUseALF()) {
    cs.picture->m_finalBuffer.createFromBuf(m_tmpBuf);
  } else if (cs.sps->getUseWrapAround())
    cs.picture->m_finalBuffer.createFromBuf(cs.picture->m_bufs[PIC_RECON_WRAP]);
  else
    cs.picture->m_finalBuffer.createFromBuf(cs.picture->m_bufs[PIC_RECONSTRUCTION]);
  if (m_threadPool) {
    initBitmap(lines * 2);
    m_rowComplete.set(0);
  }
  bool bwaitRererencePictureReady = true;
  bool bLastSlice = false;
  if (m_currDecompPic->slices.size() == 1) {
    bwaitRererencePictureReady = false;
    // single slice mode
    auto& slice = m_currDecompPic->slices[0];
    auto& bitstream = slice->parseTaskParams.bitstream;
    slice->getPic()->startProcessingTimer();
    // Decode a picture
    ITT_TASKSTART(itt_domain_prs, itt_handle_parse);
    // first slice or the previous slice not wrong
    if (slice->getCtuAddrInSlice(0) == 0) {
      m_cSliceDecoder.parseSlice(slice, &bitstream, 0, this);
      const unsigned lastCtuInSlice = slice->getCtuAddrInSlice(slice->getNumCtuInSlice() - 1);
      bLastSlice = lastCtuInSlice == slice->getPPS()->pcv->sizeInCtus - 1;
    }
    ITT_TASKEND(itt_domain_prs, itt_handle_parse);

    slice->getPic()->stopProcessingTimer();

    bitstream.clearFifo();
    bitstream.clearEmulationPreventionByteLocation();
  } else {
    for (auto& slice : m_currDecompPic->slices) {
      // parse slice task here
      auto& bitstream = slice->parseTaskParams.bitstream;
      slice->getPic()->startProcessingTimer();
      // Decode a picture
      ITT_TASKSTART(itt_domain_prs, itt_handle_parse);
      // first slice or the previous slice not wrong
      if (slice->getCtuAddrInSlice(0) == 0) {
        m_cSliceDecoder.parseSlice(slice, &bitstream, 0, nullptr);
        const unsigned lastCtuInSlice = slice->getCtuAddrInSlice(slice->getNumCtuInSlice() - 1);
        bLastSlice = lastCtuInSlice == slice->getPPS()->pcv->sizeInCtus - 1;
      }
      ITT_TASKEND(itt_domain_prs, itt_handle_parse);

      slice->getPic()->stopProcessingTimer();

      bitstream.clearFifo();
      bitstream.clearEmulationPreventionByteLocation();
    }
  }

  if (bLastSlice) {
    if (m_threadPool) {
      decodeFrameWithThreadPool();
    } else {
      // decode on frame thread
      for (int y = 0; y < lines; y++) {
        if (bwaitRererencePictureReady) waitRererencePictureReady(y);  // wait for all reference picture ready
        for (int x = 0; x < w; x++) decodeCTU(x, y, 0);
      }
    }

    if (cs.sps->getUseALF()) {
      cs.picture->m_bufs[PIC_RECONSTRUCTION].swap(m_tmpBuf);
    }

    m_currDecompPic->reconstructed = true;
  } else {
    m_currDecompPic->reconstructed = false;
  }
  m_currDecompPic->inProgress = false;
  // m_currDecompPic->done.unlock();
  m_currDecompPic->neededForOutput = m_currDecompPic->slices[0]->getPicHeader()->getPicOutputFlag();
  // here we recycle frame level tmp buffer
  Pel* residual = cs.getPictureResidualBuffer();
  m_decLib->m_frameData.recycleResidual(residual, cs.getResidualPictureSize());
  m_decLib->m_frameData.recycleCSMemory(cs.m_saoBlkParam, cs.m_ctuLoopFilterDataHorEdge, cs.m_ctuCuPtr,
                                        cs.pcv->sizeInCtus);
  cs.releaseCUTUCache();

  m_sDoneBarrier.trigger();
}

void DecLibRecon::threadProc() {
  // nothing here
  // printf("hello world!\n");
  while (m_bExitThread == false) {
    m_sEnableBarrier.wait();
    if (m_bExitThread) break;
    if (m_currDecompPic == nullptr) break;
    decodeFrame();
  }
}

void DecLibRecon::create(unsigned instanceId, DecLib* decLib, vvdec::VVDecThreadPool* pool, bool createFrameThreads) {
  // run constructor again to ensure all variables, especially in DecLibParser have been reset
  this->~DecLibRecon();
  new (this) DecLibRecon;

#if TRACE_ENABLE_ITT
  if (itt_domain_decInst.size() < instanceId + 1) {
    std::string name("DecLibRecon " + std::to_string(instanceId));
    itt_domain_decInst.push_back(__itt_domain_create(name.c_str()));
    itt_domain_decInst.back()->flags = 1;

    CHECK(itt_domain_decInst.back() != itt_domain_decInst[instanceId],
          "current decLibRecon ITT-Domain is not the last in vector. Instances created in the wrong order?");
  }
  m_itt_decInst = itt_domain_decInst[instanceId];
#endif
  m_decLib = decLib;
  // m_decodeThreadPool = threadPool;
  m_numDecThreads = std::max(1, pool ? pool->m_numWorkThreads : 1);
  m_cIntraPred = new IntraPrediction[m_numDecThreads + 1];
  m_cInterPred = new InterPrediction[m_numDecThreads + 1];
  m_cTrQuant = new TrQuant[m_numDecThreads + 1];
  m_cCuDecoder = new DecCu[m_numDecThreads + 1];
  m_cReshaper = new Reshape[m_numDecThreads + 1];
  if (createFrameThreads) m_cFrameThread = new std::thread(&DecLibRecon::threadProc, this);
  m_threadPool = pool;
  if (pool) {
    pool->addWorkProducer(this);
    m_pool = pool;
  }
  m_numRowCTUs = 0;
}

bool DecLibRecon::initBitmap(int numRows) {
    if (numRows != m_numRows) {
        delete[]m_rowDependencyBitmap;
        delete[]m_frameParseRowBitmap;
    }
    else {
        if (m_rowDependencyBitmap && m_frameParseRowBitmap) {
            for (int i = 0; i < m_numWords; i++) {
                m_rowDependencyBitmap[i] = 0;
                m_frameParseRowBitmap[i] = 0;
            }
        }
        return m_rowDependencyBitmap && m_frameParseRowBitmap;
    }
    m_numRows = numRows;

    m_numWords = (numRows + 31) >> 5;
    m_rowDependencyBitmap = new std::atomic_uint32_t[m_numWords];
    m_frameParseRowBitmap = new std::atomic_uint32_t[m_numWords];
    
    if (m_rowDependencyBitmap && m_frameParseRowBitmap) {
        for (int i = 0; i < m_numWords; i++) {
            m_rowDependencyBitmap[i] = 0;
            m_frameParseRowBitmap[i] = 0;
        }
    }
    return m_rowDependencyBitmap && m_frameParseRowBitmap;
}

void DecLibRecon::enableRow(int row) {
    uint32_t bit = 1 << (row & 31);
    std::atomic_fetch_or(&m_rowDependencyBitmap[row >> 5], bit);
    if (m_pool) wakeUpFreeThread();
}

void DecLibRecon::markRowCanBeProcessed(int row) {
    uint32_t bit = 1 << (row & 31);
    std::atomic_fetch_or(&m_frameParseRowBitmap[row >> 5], bit);
}

void DecLibRecon::doWork(int threadId) {
    unsigned long id;
    for (int w = 0; w < m_numWords; w++) {
        uint32_t oldval = m_rowDependencyBitmap[w] & m_frameParseRowBitmap[w];
        while (oldval) {
            find_first_set_bit(id, oldval);

            uint32_t bit = 1 << id;
            if (std::atomic_fetch_and(&m_rowDependencyBitmap[w], ~bit) & bit) {
                int row = w * 32 + id;
                if (row < m_numRowCTUs)
                  processRowRecon(row, threadId);
                else
                  processRowFilter(row - m_numRowCTUs, threadId);
                m_bNeedDoWork = true;
                return;
            }

            oldval = m_rowDependencyBitmap[w] & m_frameParseRowBitmap[w];
        }
    }

    m_bNeedDoWork = false;
}

void DecLibRecon::destroy() {
  if (m_cFrameThread != nullptr) {
    m_bExitThread = true;
    m_sEnableBarrier.trigger();
    m_cFrameThread->join();
    delete m_cFrameThread;
    m_cFrameThread = nullptr;
  }
  if (m_rowSync != nullptr) delete[] m_rowSync;
  m_rowSync = nullptr;
  if (m_filterRowSync != nullptr) delete[] m_filterRowSync;
  m_filterRowSync = nullptr;
  // m_decodeThreadPool = nullptr;
  m_tmpBuf.destroy();
  delete[] m_cIntraPred;
  m_cIntraPred = nullptr;
  delete[] m_cInterPred;
  m_cInterPred = nullptr;
  delete[] m_cTrQuant;
  m_cTrQuant = nullptr;
  delete[] m_cCuDecoder;
  m_cCuDecoder = nullptr;
  delete[] m_cReshaper;
  m_cReshaper = nullptr;
  free((void*)m_rowDependencyBitmap);
  free((void*)m_frameParseRowBitmap);
}

void DecLibRecon::updateFuncPtr(Picture* pcPic) {
  int bitDepth = pcPic->cs->sps->getBitDepth(CHANNEL_TYPE_LUMA);

  if (m_currFuncBitDepth == bitDepth) return;

  m_currFuncBitDepth = bitDepth;

  // init Pel Buf Ptr
  g_pelBufOP.initPelBufOps(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_BUFFER
  g_pelBufOP.initPelBufOpsX86(m_currFuncBitDepth);
#endif

  // init TransForm Ptr
  g_tCoeffOps.initTCoeffOps(m_currFuncBitDepth);
#if ENABLE_SIMD_TCOEFF_OPS
  g_tCoeffOps.initTCoeffOpsX86(m_currFuncBitDepth);
#endif

  // init SAO Ptr
  m_cSAO.initSampleAdaptiveOffset(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_SAO
  m_cSAO.initSampleAdaptiveOffsetX86(m_currFuncBitDepth);
#endif

  // init RdCost Ptr
  m_cRdCost.initRdCost(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_DIST
  m_cRdCost.initRdCostX86(m_currFuncBitDepth);
#endif

  // init LoopFilter Ptr
  m_cLoopFilter.initLoopFilter(m_currFuncBitDepth);
#if defined(TARGET_SIMD_X86) && ENABLE_SIMD_DBLF
  m_cLoopFilter.initLoopFilterX86(m_currFuncBitDepth);
#endif

  // init ALF Ptr
  m_cALF.initAdaptiveLoopFilter(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_ALF
  m_cALF.initAdaptiveLoopFilterX86(m_currFuncBitDepth);
#endif

  for (int i = 0; i < m_numDecThreads; i++) {
    // init Intra Ptr
    m_cIntraPred[i].initIntraPrediction(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_INTRAPRED
    m_cIntraPred[i].initIntraPredictionX86(m_currFuncBitDepth);
#endif

    m_cIntraPred[i].m_matrixIntraPred.initMatrixIntraPrediction(m_currFuncBitDepth);

    // init Inter Ptr
    m_cInterPred[i].initInterPrediction(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_BIO
    m_cInterPred[i].initInterPredictionX86(m_currFuncBitDepth);
#endif

    // init Quant Ptr
    m_cTrQuant[i].initQuant(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_QUANT
    m_cTrQuant[i].initQuantX86(m_currFuncBitDepth);
#endif

    // init Interpolation Ptr
    m_cInterPred[i].m_if.initInterpolationFilter(m_currFuncBitDepth);
#if ENABLE_SIMD_OPT_MCIF
    m_cInterPred[i].m_if.initInterpolationFilterX86(m_currFuncBitDepth);
#endif
  }
}

void DecLibRecon::decompressPicture(Picture* pcPic) {
  CodingStructure& cs = *pcPic->cs;

  pcPic->inProgress = true;

#ifdef TRACE_ENABLE_ITT
  // mark start of frame
  pcPic->m_itt_decLibInst = m_itt_decInst;
  __itt_frame_begin_v3(pcPic->m_itt_decLibInst, nullptr);
#endif

  // Initialise the various objects for the new set of settings
  const SPS* sps = cs.sps.get();
  const PPS* pps = cs.pps.get();

  if (m_threadPool != nullptr) {
    if (m_numRowCTUs != cs.pcv->heightInCtus) {
      if (m_rowSync != nullptr) delete[] m_rowSync;
      m_rowSync = new CTURowSync[cs.pcv->heightInCtus];
      if (m_filterRowSync != nullptr) delete[] m_filterRowSync;
      m_filterRowSync = new CTURowSync[cs.pcv->heightInCtus];
    }
    m_numRowCTUs = cs.pcv->heightInCtus;
    for (int i = 0; i < m_numRowCTUs; i++) m_rowSync[i].init();
    for (int i = 0; i < m_numRowCTUs; i++) m_filterRowSync[i].init();
  }

  for (int i = 0; i < m_numDecThreads + 1; i++) {
    if (sps->getUseReshaper()) {
      m_cReshaper[i].createDec(sps->getBitDepth(CHANNEL_TYPE_LUMA));
      m_cReshaper[i].initSlice(pcPic->slices[0]);
    }

    m_cIntraPred[i].init(sps->getChromaFormatIdc(), sps->getBitDepth(CHANNEL_TYPE_LUMA));
    m_cInterPred[i].init(&m_cRdCost, sps->getChromaFormatIdc(), sps->getMaxCUHeight(),
                         sps->getIBCFlag() ? true : false);

    // Recursive structure
    m_cTrQuant[i].init(pcPic->slices[0]);
    m_cCuDecoder[i].init(&m_cIntraPred[i], &m_cInterPred[i], &m_cReshaper[i], &m_cTrQuant[i], &m_cLoopFilter);
  }

  const int maxDepth = getLog2(sps->getMaxCUWidth()) - pps->pcv->minCUWidthLog2;
  m_cSAO.create(pps->getPicWidthInLumaSamples(), pps->getPicHeightInLumaSamples(), sps->getChromaFormatIdc(),
                sps->getMaxCUWidth(), sps->getMaxCUHeight(), maxDepth);

  if (sps->getUseALF()) {
    m_cALF.create(cs.picHeader, sps, pps, m_numDecThreads);
  }

  if (sps->getIBCFlag()) {
    const int heightInCtus = cs.pcv->heightInCtus;
    cs.initVIbcBuf(heightInCtus, sps->getChromaFormatIdc(), sps->getMaxCUHeight());
  }
  const bool doALF = cs.sps->getUseALF();
  if (doALF) {
    AdaptiveLoopFilter::preparePic(cs, &m_tmpBuf);
  }
  m_currDecompPic = pcPic;
  m_sEnableBarrier.trigger();
  if (m_cFrameThread == nullptr) {
    decodeFrame();
  }
}

Picture* DecLibRecon::waitForPrevDecompressedPic() {
  if (!m_currDecompPic) return nullptr;
  m_sDoneBarrier.wait();
  Picture* pic = m_currDecompPic;
  m_currDecompPic = nullptr;
  return pic;
}

void DecLibRecon::saoCTU(CodingStructure& cs, int x, int y, int tid) {
  const UnitArea ctuArea = getCtuArea(cs, x, y, true);
  if (cs.sps->getUseSAO()) {
    m_cSAO.SAOProcessCTU(cs, ctuArea);
  }
  if (cs.sps->getUseALF()) {
    AdaptiveLoopFilter::prepareCTU(cs, x, y);
    // apply alf here
    if (y > 0) {
      if (x > 0) {
        m_cALF.processCTU(cs, x - 1, y - 1, tid);
        cs.picture->m_rowReconCounter[y - 1].incr();
      }
      if (x == cs.pcv->widthInCtus - 1) {
        m_cALF.processCTU(cs, x, y - 1, tid);
        cs.picture->extendRowBorder(y - 1);
        cs.picture->m_rowReconCounter[y - 1].incr();
        cs.picture->m_rowCompleteCount.incr();
      }
    }
    if (y == cs.pcv->heightInCtus - 1) {
      if (x > 1) {
        m_cALF.processCTU(cs, x - 2, y, tid);
        cs.picture->m_rowReconCounter[y].incr();
      }
      if (x == cs.pcv->widthInCtus - 1) {
        if (x > 0) {
          m_cALF.processCTU(cs, x - 1, y, tid);
          cs.picture->m_rowReconCounter[y].incr();
        }
        m_cALF.processCTU(cs, x, y, tid);
        cs.picture->extendRowBorder(y);
        cs.picture->m_rowReconCounter[y].incr();
        cs.picture->m_rowCompleteCount.incr();
      }
    }
  } else {
    if (x == cs.pcv->widthInCtus - 1) {
      cs.picture->extendRowBorder(y);
      cs.picture->m_rowCompleteCount.incr();
    }
    cs.picture->m_rowReconCounter[y].incr();
  }
}

void DecLibRecon::dbfHCTU(CodingStructure& cs, int x, int y, int tid) {
  // fall back to sao
  if (y > 0) {
    if (x > 0) saoCTU(cs, x - 1, y - 1, tid);
    if (x == cs.pcv->widthInCtus - 1) saoCTU(cs, x, y - 1, tid);
  }
  if (y == cs.pcv->heightInCtus - 1) {
    if (x > 1) saoCTU(cs, x - 2, y, tid);
    if (x == cs.pcv->widthInCtus - 1) {
      if (x > 0) saoCTU(cs, x - 1, y, tid);
      saoCTU(cs, x, y, tid);
    }
  }
}

void DecLibRecon::reshapeCTU(CodingStructure& cs, int x, int y, int tid) {
  const UnitArea ctuArea = getCtuArea(cs, x, y, true);

  CtuLoopFilterData ctuLoopFilterDataVerEdge;
  memset(ctuLoopFilterDataVerEdge.lfParam, 0, sizeof(ctuLoopFilterDataVerEdge.lfParam));
  m_cLoopFilter.calcFilterStrengthsVerEdgeCTU(cs, ctuArea, ctuLoopFilterDataVerEdge.lfParam);

  if (cs.sps->getSPSTemporalMVPEnabledFlag() &&
      (cs.sps->getMaxTLayers() == 1 || cs.picture->getTLayer() < cs.sps->getMaxTLayers() - 1) &&
      !(cs.getCU(ctuArea.lumaPos(), CH_L)->slice->isIntra())) {
    m_cCuDecoder[tid].TaskDeriveDMVRMotionInfo(cs, ctuArea);
  }
  cs.picture->m_rowMotionInfoCounter[y].incr();
  m_cReshaper[tid].rspCtu(cs, x, y, 0);
  // fall back to next step
  m_cLoopFilter.loopFilterCTU(cs, MAX_NUM_CHANNEL_TYPE, x, y, EDGE_VER, ctuLoopFilterDataVerEdge.lfParam);
  if (x > 0) {
    m_cLoopFilter.loopFilterCTU(cs, MAX_NUM_CHANNEL_TYPE, x - 1, y, EDGE_HOR, nullptr);
    dbfHCTU(cs, x - 1, y, tid);
  }
  if (x == cs.pcv->widthInCtus - 1) {
    m_cLoopFilter.loopFilterCTU(cs, MAX_NUM_CHANNEL_TYPE, x, y, EDGE_HOR, nullptr);
    dbfHCTU(cs, x, y, tid);
  }
}

void DecLibRecon::waitRererencePictureReady(int row) {
  auto& cs = *m_currDecompPic->cs;
  int w = cs.pcv->widthInCtus;
  for (int x = 0; x < w; x++) {
    const UnitArea ctuArea = getCtuArea(cs, (unsigned)x, (unsigned)row, true);
    m_cCuDecoder[0].TaskWaitReferenceReady(cs, ctuArea, m_sPerLineMiHist[row],
                                           m_decLib->getNumFrameDecoder() > 1 ? true : false);
  }
}

void DecLibRecon::deriveMV(int x, int y, int tid) {
  auto& cs = *m_currDecompPic->cs;
  int w = cs.pcv->widthInCtus;
  const UnitArea ctuArea = getCtuArea(cs, (unsigned)x, (unsigned)y, true);
  m_cCuDecoder[tid].TaskWaitReferenceReady(cs, ctuArea, m_sPerLineMiHist[y],
                                           m_decLib->getNumFrameDecoder() > 1 ? true : false);
  if (x == w - 1 && m_threadPool) {
    markRowCanBeProcessed(y);
    markRowCanBeProcessed(y + m_numRowCTUs);
    if (!y) {
      m_rowSync[0].active = true;
      enableRow(y);
    }
    wakeUpFreeThread();
  }
}

void DecLibRecon::decodeCTU(int x, int y, int tid) {
  // assume all reference is available
  auto& cs = *m_currDecompPic->cs;
  const UnitArea ctuArea = getCtuArea(cs, (unsigned)x, (unsigned)y, true);
  CtuLoopFilterData& ctuLoopFilterDataHorEdge = cs.getCtuLoopFilterDataHorEdge(x, y);
  memset(ctuLoopFilterDataHorEdge.lfParam, 0, sizeof(ctuLoopFilterDataHorEdge.lfParam));
  m_cCuDecoder[tid].TaskReconAll2(cs, ctuArea);
  if (y > 0) {
    if (x > 0) reshapeCTU(cs, x - 1, y - 1, tid);
    if (x == cs.pcv->widthInCtus - 1) reshapeCTU(cs, x, y - 1, tid);
  }
  if (y == cs.pcv->heightInCtus - 1) {
    if (x > 1) reshapeCTU(cs, x - 2, y, tid);
    if (x == cs.pcv->widthInCtus - 1) {
      if (x > 0) reshapeCTU(cs, x - 1, y, tid);
      reshapeCTU(cs, x, y, tid);
    }
  }
}

bool DecLibRecon::checkRowShouldContinue(int row) {
  int numCol = m_currDecompPic->cs->pcv->widthInCtus;
  const CTURowSync& currRow = m_rowSync[row];
  int col = currRow.completed;
  if (currRow.completed == numCol) return false;
  if (row > 0 && m_rowSync[row - 1].completed != numCol && m_rowSync[row - 1].completed < col + 2) return false;
  return true;
}

void DecLibRecon::activeRow(int row) {
  CTURowSync& currRow = m_rowSync[row];
  std::lock_guard<std::mutex> guad(currRow.lock);
  if (currRow.active) return;
  if (checkRowShouldContinue(row)) {
    currRow.active = true;
    enableRow(row);
  }
}

bool DecLibRecon::checkFilterRowShouldContinue(int row) {
  int numCol = m_currDecompPic->cs->pcv->widthInCtus;
  const CTURowSync& currRow = m_filterRowSync[row];
  int col = currRow.completed;
  if (currRow.completed == numCol) return false;
  if (row > 0 && m_filterRowSync[row - 1].completed != numCol && m_filterRowSync[row - 1].completed < col + 2)
    return false;
  if (row < m_numRowCTUs - 1 && m_rowSync[row + 1].completed != numCol && m_rowSync[row + 1].completed < col + 2)
    return false;
  return true;
}

void DecLibRecon::activeFilterRow(int row) {
  CTURowSync& currRow = m_filterRowSync[row];
  std::lock_guard<std::mutex> guad(currRow.lock);
  if (currRow.active) return;
  if (checkFilterRowShouldContinue(row)) {
    currRow.active = true;
    enableRow(row + m_numRowCTUs);
  }
}

void DecLibRecon::processRowFilter(int row, int threadId) {
  CTURowSync& r = m_filterRowSync[row];
  {
    std::lock_guard<std::mutex> guad(r.lock);
    if (!r.active) return;
    if (r.busy) return;
    r.busy = true;
  }
  auto& cs = *m_currDecompPic->cs;
  int numCol = cs.pcv->widthInCtus;
  int numRow = cs.pcv->heightInCtus;
  while (r.completed < numCol) {
    reshapeCTU(cs, r.completed, row, threadId);
    r.completed++;
    if (r.completed >= 2 && row + 1 < numRow) {
      activeFilterRow(row + 1);
    }
    if (r.completed < numCol - 1) {
      std::lock_guard<std::mutex> guad(r.lock);
      if (checkFilterRowShouldContinue(row) == false) {
        r.busy = false;
        r.active = false;
        return;
      }
    }
  }
  r.busy = false;
  r.active = false;
  m_rowComplete.incr();
  if (row == numRow - 1) {
    int v = m_rowComplete.get();
    while (v != numRow) {
      v = m_rowComplete.waitForChange(v);
    }
    m_allRowDone.trigger();
  }
}

void DecLibRecon::processRowRecon(int row, int threadId) {
  CTURowSync& r = m_rowSync[row];
  {
    std::lock_guard<std::mutex> guad(r.lock);
    if (!r.active) return;
    if (r.busy) return;
    r.busy = true;
  }
  auto& cs = *m_currDecompPic->cs;
  int numCol = cs.pcv->widthInCtus;
  int numRow = cs.pcv->heightInCtus;
  while (r.completed < numCol) {
    // decodeCTU(r.completed, row, threadId);
    {
      const UnitArea ctuArea = getCtuArea(cs, (unsigned)r.completed, (unsigned)row, true);
      CtuLoopFilterData& ctuLoopFilterDataHorEdge = cs.getCtuLoopFilterDataHorEdge(r.completed, row);
      memset(ctuLoopFilterDataHorEdge.lfParam, 0, sizeof(ctuLoopFilterDataHorEdge.lfParam));
      m_cCuDecoder[threadId].TaskReconAll2(cs, ctuArea);
    }
    // active next row
    r.completed++;
    if (r.completed >= 2 && row + 1 < numRow) {
      activeRow(row + 1);
    }
    if (row > 0 && r.completed >= 2) {
      activeFilterRow(row - 1);
    }
    // detect we should wait for last row here
    if (r.completed < numCol - 1) {
      std::lock_guard<std::mutex> guad(r.lock);
      if (checkRowShouldContinue(row) == false) {
        r.busy = false;
        r.active = false;
        return;
      }
    }
  }
  r.busy = false;
  r.active = false;
}

void DecLibRecon::decodeFrameWithThreadPool() {
  bool bEnableRowEarly = m_currDecompPic->slices.size() == 1;
  if (bEnableRowEarly == false) {
    int lines = m_currDecompPic->cs->pcv->heightInCtus;
    for (int y = 0; y < lines; y++) {
      waitRererencePictureReady(y);
      markRowCanBeProcessed(y);
      markRowCanBeProcessed(y + m_numRowCTUs);
      if (!y) {
        m_rowSync[0].active = true;
        enableRow(y);
      }
      wakeUpFreeThread();
    }
  }
  static const int block_ms = 250;
  while (m_allRowDone.timedWait(block_ms)) {
    wakeUpFreeThread();
  }
}
//! \}
