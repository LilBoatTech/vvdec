#ifndef VVDEC_THREADPOOL_H
#define VVDEC_THREADPOOL_H
#include <mutex>
#include <vector>
#include "threading.h"

namespace vvdec {

class VVDecThreadPool;
class VVDecThread;

#ifdef __GNUC__
#  define find_first_set_bit(id, x) id = (unsigned long)__builtin_ctz(x)
#elif defined(_MSC_VER)
#  include <intrin.h>
#  define find_first_set_bit(id, x) _BitScanForward(&id, x)

#endif

typedef uint32_t ThreadBitmap;
constexpr int MAX_POOL_THREADS = 32;

class VVDecWorkProducer {
 public:
  VVDecThreadPool* m_pool;
  bool m_bNeedDoWork;
  VVDecWorkProducer() : m_pool(NULL), m_bNeedDoWork(false) {}
  virtual ~VVDecWorkProducer() {}
  virtual void doWork(int workerThreadId) = 0;
  void wakeUpFreeThread();
};

class VVDecThreadPool {
 public:
  std::atomic_uint32_t m_freeThreads;
  int m_numWorkThreads;
  bool m_activeFlag;
  std::vector<VVDecWorkProducer*> m_wps;
  VVDecThread* m_workThreads;

  VVDecThreadPool();
  ~VVDecThreadPool();
  // should be called before the threads pool start
  void addWorkProducer(VVDecWorkProducer* pWorkProducer) {
	  m_wps.push_back(pWorkProducer);
  }
  bool create(int numThreads, int maxProviders);
  bool start();
  void stopWorkers();
  static VVDecThreadPool* allocThreadPools(int maxProvider, int numPoolThreads, int& numPool);
};
}  // namespace vvdec

#endif  // ifndef VVDEC_THREADPOOL_H
