/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

/* $Id: random_choice.h 693 2004-03-16 15:48:04Z wistaria $ */

#ifndef LOOPER_RANDOM_CHOICE_H
#define LOOPER_RANDOM_CHOICE_H

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>
#include <boost/config.hpp>
#include <boost/random/uniform_01.hpp>

// Define the following macro if you want the original initialization
// routine of O(N^2).
// #define USE_INIT_WALKER1977

namespace looper {

template<class IntType = int, class RealType = double>
class random_choice
{
public:
  typedef RealType input_type;
  typedef IntType result_type;

  random_choice() : n_(0) {}
  template<class CONT>
  explicit random_choice(const CONT& weights)
  {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif

#ifndef USE_INIT_WALKER1977
    init(weights);
#else
    init_walker1977(weights);
#endif
  }
  // compiler-generated copy ctor and assignment operator are fine

  template<class Engine>
  result_type operator()(Engine& eng) const
  {
    result_type x = result_type(RealType(n_) * eng());
    return (eng() < cutoff(x)) ? x : alias(x);
  }

  // Initialization routine with complexity O(N).
  template<class CONT>
  void init(const CONT& weights)
  {
    assert(weights.size() != 0);

    n_ = weights.size();
    RealType norm = RealType(0);
    for (result_type i = 0; i < n_; ++i) norm += weights[i];
    norm = n_ / norm;

    // Initialize arrays.  We will reorder the elements in `array', so
    // that all the negative elements precede the positive ones.
    table_.resize(n_);
    std::vector<std::pair<RealType, result_type> > array(n_);
    typename std::vector<std::pair<RealType, result_type> >::iterator
      neg_p = array.begin();
    typename std::vector<std::pair<RealType, result_type> >::iterator
      pos_p = array.end();
    for (result_type i = 0; i < n_; ++i) {
      RealType b = norm * weights[i] - RealType(1);
      if (b < RealType(0)) {
        *neg_p = std::make_pair(b, i);
        ++neg_p;
      } else {
        --pos_p;
        *pos_p = std::make_pair(b, i);
      }
    }

    // Note: now `pos_p' is pointing the first non-negative element in
    // the array.

    // Assign alias and cutoff values
    for (neg_p = array.begin(); neg_p != array.end(); ++neg_p) {
      cutoff(neg_p->second) = 1 + neg_p->first;
      alias(neg_p->second) = pos_p->second;
      pos_p->first += neg_p->first;
      if (pos_p->first < 0) ++pos_p;
    }
  }

  // Original O(N^2) initialization routine given in A. W. Walker, ACM
  // Trans. Math. Software, 3, 253 (1077).
  template<class CONT>
  void init_walker1977(const CONT& weights, RealType tol = 1.0e-10) {
    assert(weights.size() > 0);

    n_ = weights.size();
    RealType norm = 0;
    for (result_type i = 0; i < n_; ++i) norm += weights[i];
    norm = n_ / norm;

    // Initialize arrays.
    table_.assign(n_, std::make_pair(1, 0));
    std::vector<RealType> b(n_);
    for (result_type i = 0; i < n_; ++i) {
      cutoff(i) = 1;
      alias(i) = i;
      b[i] = norm * weights[i] - 1;
    }

    for (result_type i = 0; i < n_; ++i) {
      // Find the largest positive and negative differences and their
      // positions.
      RealType sum = 0;
      RealType minval = 0;
      RealType maxval = 0;
      result_type minpos, maxpos;
      for (result_type j = 0; j < n_; ++j) {
        sum += std::abs(b[j]);
        if (b[j] <= minval) {
          minval = b[j];
          minpos = j;
        }
        if (b[j] >= maxval) {
          maxval = b[j];
          maxpos = j;
        }
      }

      if (sum < tol) break;

      // Assign alias and cutoff values
      cutoff(minpos) = 1 + minval;
      alias(minpos) = maxpos;
      b[maxpos] += minval;
      b[minpos] = 0;
    }
  }

  template<class CONT>
  bool check(const CONT& weights, RealType tol = 1.0e-10) const {
    bool r = true;
    tol *= n_;

    RealType norm = 0;
    for (result_type i = 0; i < n_; ++i) norm += weights[i];
    norm = n_ / norm;

    for (result_type i = 0; i < n_; ++i) {
      RealType p = cutoff(i);
      for (result_type j = 0; j < n_; ++j)
        if (alias(j) == i) p += 1 - cutoff(j);
      if (std::abs(p - norm * weights[i]) > tol) r = false;
    }
    return r;
  }

protected:
  RealType& cutoff(result_type i) { return table_[i].first; }
  const RealType& cutoff(result_type i) const { return table_[i].first; }

  result_type& alias(result_type i) { return table_[i].second; }
  const result_type& alias(result_type i) const { return table_[i].second; }

private:
  result_type n_; // number of choices
  std::vector<std::pair<RealType, result_type> > table_;
                  // first element:  cutoff value
                  // second element: alias
};

} // end namespace looper

#endif // LOOPER_RANDOM_CHOICE_H
