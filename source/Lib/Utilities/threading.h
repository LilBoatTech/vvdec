#ifndef VVDEC_THREADING_H
#define VVDEC_THREADING_H

#include <chrono>
#include <condition_variable>
#include <mutex>
namespace vvdec {
class VVDecCondVar {
 public:
  VVDecCondVar() { m_counter = 0; }

  ~VVDecCondVar() {}

  void wait() {
    std::unique_lock<std::mutex> lk(cs);
    while (m_counter == 0) cv.wait(lk);
    m_counter--;
  }

  bool timedWait(uint32_t milliseconds) {
    bool bTimeOut = false;
    std::unique_lock<std::mutex> lk(cs);
    while (m_counter == 0 && bTimeOut == false) {
      bTimeOut = cv.wait_for(lk, std::chrono::milliseconds(milliseconds)) == std::cv_status::timeout;
    }
    if (m_counter > 0) {
      bTimeOut = false;
      m_counter--;
    }
    return bTimeOut;
  }

  void trigger() {
    {
      std::unique_lock<std::mutex> lk(cs);
      if (m_counter < 0xffffffff) m_counter++;
    }
    cv.notify_all();
  }

 protected:
  std::condition_variable cv;
  std::mutex cs;
  uint32_t m_counter;
};

class VVDecThreadCounter {
 public:
  VVDecThreadCounter() { m_val = 0; }
  ~VVDecThreadCounter() {}
  int waitForChange(int prev) {
    std::unique_lock<std::mutex> lk(m_cs);
    while (m_val == prev) m_cv.wait(lk);
    return m_val;
  }
  int get() {
    std::unique_lock<std::mutex> lk(m_cs);
    int ret = m_val;
    return ret;
  }
  void set(int newval) {
    {
      std::unique_lock<std::mutex> lk(m_cs);
      m_val = newval;
    }
    m_cv.notify_all();
  }
  void incr() {
      add(1);
  }
  void add(int d) {
    {
      std::unique_lock<std::mutex> lk(m_cs);
      m_val += d;
    }
    m_cv.notify_all();
  }

 protected:
  std::mutex m_cs;
  std::condition_variable m_cv;
  int m_val;
};
}  // namespace vvdec

#endif  // ifndef VVDEC_THREADING_H
