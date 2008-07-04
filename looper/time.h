/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TIME_H
#define LOOPER_TIME_H

#include <alps/config.h>
#include <alps/osiris.h>
#include <cmath>
#include <complex>

namespace looper {

inline std::complex<double> ctime(double t) {
  return std::exp(std::complex<double>(0, 2*M_PI*t));
}

template<bool STORE_CTIME = true>
class imaginary_time;

// default version
// store complex time inside the object
template<bool STORE_CTIME>
class imaginary_time {
public:
  explicit imaginary_time(double t) : time_(t), ctime_(::looper::ctime(t)) {}
  operator double() const { return time_; }
  std::complex<double> const& ctime() const { return ctime_; }
private:
  double time_;
  std::complex<double> ctime_;
};

// another version
// calculate complex time when needed
template<>
class imaginary_time<false> {
public:
  explicit imaginary_time(double t) : time_(t) {}
  operator double() const { return time_; }
  std::complex<double> ctime() const { return ::looper::ctime(time_); }
private:
  double time_;
};

inline std::complex<double> const& ctime(imaginary_time<true> const& t) {
  return t.ctime();
}

inline std::complex<double> ctime(imaginary_time<false> const& t) {
  return t.ctime();
}

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<bool F>
alps::ODump& operator<<(alps::ODump& dp, looper::imaginary_time<F> const& t) {
  dp << static_cast<double>(t); return dp;
}

template<bool F>
alps::IDump& operator>>(alps::IDump& dp, looper::imaginary_time<F>& t) {
  double t_in;
  dp >> t_in;
  t = looper::imaginary_time<F>(t_in);
  return dp;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_TIME_H
