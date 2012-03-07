/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_CAPACITY_H
#define LOOPER_CAPACITY_H

#if defined(LOOPER_ENABLE_OPENMP) && defined(_OPENMP) && !defined(LOOPER_OPENMP)
# define LOOPER_OPENMP
#endif

#include <vector>
#include <iostream>
#ifdef LOOPER_OPENMP
# include <omp.h>
#endif

namespace looper {

class vector_capacity {
public:
#ifndef LOOPER_OPENMP
  template<typename OP, typename EST, typename FRAG>
  vector_capacity(std::vector<double> const& times, std::vector<OP> const& operators,
    std::vector<OP> const& operators_p, std::vector<EST> const& estimates,
    std::vector<FRAG> const& fragments) {
    capacity_[0] = times.capacity();
    capacity_[1] = std::max(operators.capacity(), operators_p.capacity());
    capacity_[2] = estimates.capacity();
    capacity_[3] = fragments.capacity();
  }
  template<typename OP, typename ESTI, typename ESTN, typename FRAG>
  vector_capacity(std::vector<double> const& times, std::vector<OP> const& operators,
    std::vector<OP> const& operators_p, std::vector<ESTI> const& estimates_i,
    std::vector<ESTN> const& estimates_n, std::vector<FRAG> const& fragments) {
    capacity_[0] = times.capacity();
    capacity_[1] = std::max(operators.capacity(), operators_p.capacity());
    capacity_[2] = std::max(estimates_i.capacity(), estimates_n.capacity());
    capacity_[3] = fragments.capacity();
  }
#else
  template<typename OP, typename EST, typename FRAG>
  vector_capacity(std::vector<std::vector<double> > const& times_g,
    std::vector<std::vector<OP> > const& operators_g,
    std::vector<std::vector<OP> > const& operators_pg,
    std::vector<std::vector<EST> > const& estimates_g,
    std::vector<FRAG> const& fragments) {
    for (int i = 0; i < 4; ++i) capacity_[i] = 0;
    for (int tid = 0; tid < omp_get_max_threads(); ++tid) {
      capacity_[0] = std::max(capacity_[0], times_g[tid].capacity());
      capacity_[1] = std::max(capacity_[1], operators_g[tid].capacity());
      capacity_[1] = std::max(capacity_[1], operators_pg[tid].capacity());
      capacity_[2] = std::max(capacity_[2], estimates_g[tid].capacity());
    }
    capacity_[3] = fragments.capacity();
  }
  template<typename OP, typename ESTI, typename ESTN, typename FRAG>
  vector_capacity(std::vector<std::vector<double> > const& times_g,
    std::vector<std::vector<OP> > const& operators_g,
    std::vector<std::vector<OP> > const& operators_pg,
    std::vector<std::vector<ESTI> > const& estimates_ig,
    std::vector<std::vector<ESTN> > const& estimates_ng,
    std::vector<FRAG> const& fragments) {
    for (int i = 0; i < 4; ++i) capacity_[i] = 0;
    for (int tid = 0; tid < omp_get_max_threads(); ++tid) {
      capacity_[0] = std::max(capacity_[0], times_g[tid].capacity());
      capacity_[1] = std::max(capacity_[1], operators_g[tid].capacity());
      capacity_[1] = std::max(capacity_[1], operators_pg[tid].capacity());
      capacity_[2] = std::max(capacity_[2], estimates_ig[tid].capacity());
      capacity_[2] = std::max(capacity_[2], estimates_ng[tid].capacity());
    }
    capacity_[3] = fragments.capacity();
  }
#endif
  void report() const {
    std::cerr << "Info: vector capacity (max): times = " << capacity_[0]
              << " operators = " << capacity_[1]
              << " estimates = " << capacity_[2]
              << " fragments = " << capacity_[3] << std::endl;
  }

private:
  std::size_t capacity_[4]; // times, operators, estimates, fragments
};

}

#endif // LOOPER_CAPACITY_H
