/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2009-2010 by Synge Todo <wistaria@comp-phys.org>,
*                            Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef LOOPER_ATOMIC_H
#define LOOPER_ATOMIC_H

//
// compare-and-swap
//

namespace looper {
bool compare_and_swap(int& variable, int oldval, int newval);
}

#if defined(__linux__) && defined(__x86_64__)

// Linux on X86_64
#if (__GNUC__ >= 4) && (__GNUC_MINOR__ >= 1) && !defined(__ICC)

// for GCC 4.1 or later
namespace looper {
inline bool compare_and_swap(int& variable, int oldval, int newval) {
  return __sync_bool_compare_and_swap(&variable, oldval, newval);
}
}

#else

// taken from http://stackoverflow.com/questions/833122/xchg-example-for-64-bit-integer
namespace looper {
inline bool compare_and_swap(int& variable, int oldval, int newval) {
  // int *ptr = &variable;
  unsigned char ret;
  __asm__ __volatile__ (
                        " lock\n"
                        " cmpxchgl %2,%1\n"
                        " sete %0\n"
                        : "=q" (ret), "=m" (variable)
                        : "r" (newval), "m" (variable), "a" (oldval)
                        : "memory");
  return ret;
}
}

#endif

#elif defined(__APPLE__)

// for Mac OS X
#include <libkern/OSAtomic.h>
namespace looper {
inline bool compare_and_swap(int& variable, int oldval, int newval) {
  return OSAtomicCompareAndSwap32(oldval, newval, &variable);
}
}

#elif defined(__FCC_VERSION)

// for Solaris
#include <atomic.h>
#include <boost/static_assert.hpp>
namespace looper {
inline bool compare_and_swap(int& variable, int oldval, int newval) {
  BOOST_STATIC_ASSERT(sizeof(int) == sizeof(uint_t));
  int res = atomic_cas_uint((uint_t*)&variable, oldval, newval);
  return res == oldval;
}
}

#else

#error "atomic operations are not defined for the present architecture"

#endif

//
// thread-safe counter
//

namespace looper {
class atomic_counter {
public:
  atomic_counter(int init = 0) : counter_(init) {}
  void reset(int init = 0) { counter_ = init; }
  int operator++(int) {
    while (true) {
      int val = counter_;
      if (compare_and_swap(counter_, val, val+1)) return val;
    }
  }
  int operator++() {
    while (true) {
      int val = counter_;
      if (compare_and_swap(counter_, val, val+1)) return val+1;
    }
  }
  int non_atomic_increment() { return counter_++; }
  int operator()() const { return counter_; }
private:
  int counter_;
};
}

#endif // LOOPER_ATOMIC_H
