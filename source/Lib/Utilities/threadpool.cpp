#if defined(_MSC_VER)
#  include <Windows.h>
#endif
#include <stdlib.h>  // for malloc & free
#include <string.h>
#include <thread>
#include <algorithm>
#include "threadpool.h"
#include "threading.h"

#include <new>
namespace vvdec {

class VVDecThread {
 private:
  VVDecThreadPool& m_pool;
  int m_threadId;
  VVDecCondVar m_wakeCondVar;
  std::thread* thread;
  VVDecThread& operator=(const VVDecThread&);

 public:
  VVDecWorkProducer* m_currWorkProducer;

  VVDecThread(VVDecThreadPool& pool, int id) : m_pool(pool), m_threadId(id) {
    thread = new std::thread(&VVDecThread::threadEntry, this);
  }
  ~VVDecThread() {}

  void threadEntry();
  void wakeup() { m_wakeCondVar.trigger(); }
  void stop() {
    if (thread) {
      ThreadBitmap threadMask = (ThreadBitmap)1 << m_threadId;
      while (!(m_pool.m_freeThreads & threadMask))
        std::this_thread::sleep_for(std::chrono::microseconds(1));
      wakeup();
      thread->join();
      delete thread;
      thread = nullptr;
    }
  }
};

void VVDecThread::threadEntry() {
  const ThreadBitmap threadMask = (ThreadBitmap)1 << m_threadId;
  std::atomic_fetch_or(&m_pool.m_freeThreads, threadMask);
  m_currWorkProducer = nullptr;
  m_wakeCondVar.wait();
  while (m_pool.m_activeFlag) {
    do {
      m_currWorkProducer->doWork(m_threadId);
      if (!m_currWorkProducer->m_bNeedDoWork) {
          auto it = std::find_if(m_pool.m_wps.begin(), m_pool.m_wps.end(),
              [](VVDecWorkProducer* job) {return job->m_bNeedDoWork; });
          if (it != m_pool.m_wps.end()) {
              m_currWorkProducer = *it;
          }
       }
    } while (m_currWorkProducer->m_bNeedDoWork);
    std::atomic_fetch_or(&m_pool.m_freeThreads, threadMask);
    m_wakeCondVar.wait();
  }
  std::atomic_fetch_or(&m_pool.m_freeThreads, threadMask);
}

void VVDecWorkProducer::wakeUpFreeThread() {
  ThreadBitmap cache = m_pool->m_freeThreads;
  while (cache) {
    unsigned long threadId;
    find_first_set_bit(threadId, cache);
    ThreadBitmap threadMask = (ThreadBitmap)1 << threadId;
    if (std::atomic_fetch_and(&m_pool->m_freeThreads, ~threadMask) & threadMask) {
      VVDecThread& worker = m_pool->m_workThreads[threadId];
      worker.m_currWorkProducer = this;
      worker.wakeup();
      return;
    }
    cache = m_pool->m_freeThreads;
  }
  m_bNeedDoWork = true;
}

VVDecThreadPool* VVDecThreadPool::allocThreadPools(int maxProvider, int numPoolThreads, int& numPools) {
  numPoolThreads = numPoolThreads <= 0 ? std::thread::hardware_concurrency() : numPoolThreads;
  numPools = (numPoolThreads + MAX_POOL_THREADS - 1) / MAX_POOL_THREADS;
  if (!numPools) numPools = 1;
  VVDecThreadPool* pools = new VVDecThreadPool[numPools];
  if (pools) {
    for (int i = 0; i < numPools; i++) {
      int numThreads = std::min<int>(MAX_POOL_THREADS, numPoolThreads);
      if (!pools[i].create(numThreads, maxProvider)) {
        delete[] pools;
        numPools = 0;
        return nullptr;
      }
      numPoolThreads -= numPoolThreads;
    }
  } else {
    numPools = 0;
  }
  return pools;
}

VVDecThreadPool::VVDecThreadPool() {}

bool VVDecThreadPool::create(int numThreads, int maxProviders) {
  m_numWorkThreads = numThreads;
  m_workThreads = (VVDecThread*)malloc(sizeof(VVDecThread) * numThreads);
  m_freeThreads = 0;
  if (m_workThreads)
    for (int i = 0; i < numThreads; i++) new (m_workThreads + i) VVDecThread(*this, i);
  m_wps.reserve(maxProviders);

  return m_workThreads ;
}

bool VVDecThreadPool::start() {
  m_activeFlag = true;
  return true;
}

void VVDecThreadPool::stopWorkers() {
  if (m_workThreads) {
    m_activeFlag = false;
    for (int i = 0; i < m_numWorkThreads; i++) {
      // wait for thread to exit & join
      m_workThreads[i].stop();
    }
  }
}

VVDecThreadPool::~VVDecThreadPool() {
  free(m_workThreads);
}
}  // namespace vvdec
