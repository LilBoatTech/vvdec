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

/** \file     DecLibRecon.h
    \brief    decoder class (header)
*/

#ifndef DECLIB_RECON_H
#define DECLIB_RECON_H
#include <thread>  // using std::thread
#include "DecLib.h"
#include "CommonLib/CommonDef.h"
#include "CommonLib/Picture.h"
#include "CommonLib/RdCost.h"
#include "CommonLib/Reshape.h"
#include "CommonLib/LoopFilter.h"
#include "CommonLib/AdaptiveLoopFilter.h"
#include "CommonLib/SampleAdaptiveOffset.h"

//#include "Utilities/NoMallocThreadPool.h"
#include "Utilities/threadpool.h"
#include "Utilities/threading.h"

class DecLibRecon;
class IntraPrediction;
class InterPrediction;
class TrQuant;
class DecCu;

//! \ingroup DecoderLib
//! \{
// ====================================================================================================================
// Class definition
// ====================================================================================================================
using vvdec::VVDecCondVar;
struct CTURowSync {
  std::mutex lock;
  volatile bool active;
  volatile bool busy;
  volatile uint32_t completed;
  void init() {
    active = false;
    busy = false;
    completed = 0;
  }
};

/// decoder class
class DecLibRecon : public vvdec::VVDecWorkProducer {
    std::atomic_uint32_t* m_rowDependencyBitmap; 
    std::atomic_uint32_t* m_frameParseRowBitmap;
    //RWLock m_mutex;

    // number of words in the bitmap
    int m_numWords;

    int m_numRows;

public:
    // called before decoding each frame
    bool initBitmap(int numRows);

   
    void enableRow(int row);
    void markRowCanBeProcessed(int row);
    void doWork(int threadId);
 private:
  // functional classes
  IntraPrediction* m_cIntraPred = nullptr;
  InterPrediction* m_cInterPred = nullptr;
  TrQuant* m_cTrQuant = nullptr;
  DecCu* m_cCuDecoder = nullptr;
  RdCost m_cRdCost;
  Reshape* m_cReshaper = nullptr;  ///< reshaper class
  std::thread* m_cFrameThread = nullptr;
  LoopFilter m_cLoopFilter;
  SampleAdaptiveOffset m_cSAO;
  AdaptiveLoopFilter m_cALF;

  int m_numDecThreads = 0;
  // NoMallocThreadPool* m_decodeThreadPool;
  vvdec::VVDecThreadPool* m_threadPool;
  PelStorage m_tmpBuf;
  DecLib* m_decLib = nullptr;

  Picture* m_currDecompPic = nullptr;
#if TRACE_ENABLE_ITT
  __itt_domain* m_itt_decInst = nullptr;
#endif

  int m_currFuncBitDepth;
  // int m_numRows;
  int m_numRowCTUs = 0;
  CTURowSync* m_rowSync = nullptr;
  CTURowSync* m_filterRowSync = nullptr;

 public:
  DecLibRecon();
  ~DecLibRecon() = default;
  DecLibRecon(const DecLibRecon&) = delete;
  DecLibRecon(const DecLibRecon&&) = delete;

  void create(unsigned instanceId, DecLib* decLib, vvdec::VVDecThreadPool* pool, bool createFrameThreads);
  void destroy();

  void updateFuncPtr(Picture* pcPic);
  void decompressPicture(Picture* pcPic);
  Picture* waitForPrevDecompressedPic();
  Picture* getCurrPic() const { return m_currDecompPic; }
  void deriveMV(int x, int y, int tid);

 private:
  void borderExtPic(Picture* pic);
  void threadProc();
  DecSlice m_cSliceDecoder;
  bool m_bExitThread = false;
  VVDecCondVar m_sEnableBarrier;
  VVDecCondVar m_sDoneBarrier;
  VVDecCondVar m_allRowDone;
  vvdec::VVDecThreadCounter m_rowComplete;
  std::vector<MotionHist> m_sPerLineMiHist;
  void waitRererencePictureReady(int row);
  void decodeCTU(int x, int y, int tid);
  void reshapeCTU(CodingStructure& cs, int x, int y, int tid);
  void dbfHCTU(CodingStructure& cs, int x, int y, int tid);
  void saoCTU(CodingStructure& cs, int x, int y, int tid);
  void processRowRecon(int row, int threadId);
  void processRowFilter(int row, int threadId);
  void activeRow(int row);
  void activeFilterRow(int row);
  bool checkRowShouldContinue(int row);
  bool checkFilterRowShouldContinue(int row);
  void decodeFrameWithThreadPool();
  void decodeFrame();
};

//! \}

#endif  // DECLIB_RECON_H
