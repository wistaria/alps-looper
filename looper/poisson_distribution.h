/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_POISSON_DISTRIBUTION_H
#define LOOPER_POISSON_DISTRIBUTION_H

#include <cmath>
#include <cassert>
#include <iostream>
#include <math.h>
#include <boost/limits.hpp>
#include <boost/static_assert.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

//
// generating random numbers according to Poisson distribution
//

template<typename IntType = int, typename RealType = double, int THRESHOLD = 16>
class poisson_distribution {
public:
  typedef RealType input_type;
  typedef IntType result_type;

  explicit poisson_distribution(const RealType& mean = RealType(1)) : mean_(mean) {
#ifndef BOOST_NO_STDC_NAMESPACE
    // allow for Koenig lookup
    using std::exp;
    using std::log;
    using std::sqrt;
#endif
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif
    if (mean_ <= RealType(0))
      boost::throw_exception(std::invalid_argument("poisson_distribution"));
    exp_mean_ = exp(-mean_);
    if (mean_ >= RealType(THRESHOLD)) {
      sqr_ = sqrt(RealType(2) * mean_);
      alxm_ = log(mean_);
      gm_ = mean_ * alxm_ - lgamma(mean_ + RealType(1));
    }
  }
  // compiler-generated copy ctor and assignment operator are fine

  RealType mean() const { return mean_; }
  void reset() {}

  template<class Engine>
  result_type operator()(Engine& eng) {
#ifndef BOOST_NO_STDC_NAMESPACE
    // allow for Koenig lookup
    using std::exp;
    using std::tan;
#endif
    if (mean_ < RealType(THRESHOLD)) {
      // Use O(mean_) method
      RealType product = RealType(1);
      for (result_type m = 0; ; ++m) {
        product *= eng();
        if (product <= exp_mean_)
          return m;
      }
    } else {
      double em;
      double y;
      double t;
      do {
        do {
          y = tan(M_PI * eng());
          em = sqr_ * y + mean_;
        } while (em < 0.0);
        em = floor(em);
        t = 0.9 * (RealType(1) + y * y) * exp(em * alxm_ - lgamma(em + RealType(1)) - gm_);
      } while (eng() > t);
      return em;
    }
  }

private:
  RealType mean_;
  // some precomputed data from the parameters
  RealType exp_mean_;
  RealType sqr_;
  RealType alxm_;
  RealType gm_;
};

} // end namespace looper

#endif // LOOPER_POISSON_DISTRIBUTION_H
