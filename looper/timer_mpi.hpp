/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2010-2012 by Synge Todo <wistaria@comp-phys.org>,
*                            Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>,
*                            Hideyuki Shitara <shitara.hide@jp.fujitsu.com>
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

#ifndef ALPS_PARAPACK_TIMER_MPI_HPP
#define ALPS_PARAPACK_TIMER_MPI_HPP

// ALPS_ENBALE_TIMER
// ALPS_ENABLE_TIMER_TRACE
// ALPS_ENABLE_TIMER_DETAILED
// ALPS_ENABLE_TIMER_BARRIER

#include "timer.hpp"
#include <mpi.h>

#ifdef ALPS_ENABLE_TIMER

namespace alps {
namespace parapack {
namespace detail {

struct clock_mpi {
  static double get_time() {
    return MPI_Wtime();
  }
};

template<class CLOCK>
class timer_mpi : public timer_base<CLOCK> {
private:
  typedef timer_base<CLOCK> base_type;
public:
  BOOST_STATIC_CONSTANT(int, start_barrier = (1 << 1));
  BOOST_STATIC_CONSTANT(int, stop_barrier = (1 << 2));
  BOOST_STATIC_CONSTANT(int, barrier = (start_barrier | stop_barrier));
  timer_mpi(MPI_Comm comm) : base_type(), comm_(comm) {}
  #ifdef ALPS_ENABLE_TIMER_BARRIER
  void registrate(std::size_t id, std::string const& label, int option = 0) {
    base_type::registrate(id, label, option);
    if (id >= barrier_.size()) barrier_.resize(id + 1);
    barrier_[id] = option;
  }
  void start(std::size_t id, bool fjtool = true) {
    if (barrier_[id] & start_barrier) MPI_Barrier(comm_);
    base_type::start(id, fjtool);
  }
  void stop(std::size_t id, bool fjtool = true) {
    if (barrier_[id] & stop_barrier) MPI_Barrier(comm_);
    base_type::stop(id, fjtool);
  }
  #endif
  void summarize(std::ostream& os = std::clog) const {
    int np, rank;
    MPI_Comm_size(comm_, &np);
    MPI_Comm_rank(comm_, &rank);
    int size = base_type::get_labels().size();
    counts_tmp_.resize(size);
    counts_tmp_ = base_type::get_counts();
    counts_ave_.resize(size);
    counts_min_.resize(size);
    counts_max_.resize(size);
    sums_tmp_.resize(size);
    sums_tmp_ = base_type::get_measurements();
    sums_ave_.resize(size);
    sums_min_.resize(size);
    sums_max_.resize(size);
#ifdef __linux
    vm_tmp_.resize(0); // Peak and Hwm
    vm_ave_.resize(2);
    vm_min_.resize(2);
    vm_max_.resize(2);
    std::map<std::string,int> vm_info_ = base_type::vmem_info();
    std::vector<std::string> vm_labels_(0);
    for (std::map<std::string,int>::iterator ivm = vm_info_.begin();
      ivm != vm_info_.end(); ++ivm) {
      vm_labels_.push_back(ivm->first);
      vm_tmp_.push_back(ivm->second);
    }
    MPI_Reduce(&vm_tmp_[0], &vm_ave_[0], 2, MPI_INT, MPI_SUM, 0, comm_);
    MPI_Reduce(&vm_tmp_[0], &vm_min_[0], 2, MPI_INT, MPI_MIN, 0, comm_);
    MPI_Reduce(&vm_tmp_[0], &vm_max_[0], 2, MPI_INT, MPI_MAX, 0, comm_);
#endif
    MPI_Reduce(&counts_tmp_[0], &counts_ave_[0], size, MPI_DOUBLE, MPI_SUM, 0, comm_);
    MPI_Reduce(&counts_tmp_[0], &counts_min_[0], size, MPI_DOUBLE, MPI_MIN, 0, comm_);
    MPI_Reduce(&counts_tmp_[0], &counts_max_[0], size, MPI_DOUBLE, MPI_MAX, 0, comm_);
    MPI_Reduce(&sums_tmp_[0], &sums_ave_[0], size, MPI_DOUBLE, MPI_SUM, 0, comm_);
    MPI_Reduce(&sums_tmp_[0], &sums_min_[0], size, MPI_DOUBLE, MPI_MIN, 0, comm_);
    MPI_Reduce(&sums_tmp_[0], &sums_max_[0], size, MPI_DOUBLE, MPI_MAX, 0, comm_);

    counts_ave_ /= np;
    sums_ave_ /= np;

#ifdef __linux
    for (int i = 0; i < vm_ave_.size(); ++i) {
      vm_ave_[i] /= np;
    }
#endif

    if (rank == 0) {
      os << "timer: enabled\n"
#ifdef ALPS_ENABLE_TIMER_TRACE
         << "timer: trace = enabled\n"
#else
         << "timer: trace = disabled\n"
#endif
#ifdef ALPS_ENABLE_TIMER_DETAILED
         << "timer: detailed report = enabled\n"
#else
         << "timer: detailed report = disabled\n"
#endif
#ifdef ALPS_ENABLE_TIMER_BARRIER
         << "timer: barrier synchronization = enabled"
#else
         << "timer: barrier synchronization = disabled"

#endif
         << std::endl;
#ifdef __linux
      for (int i = 0; i < 2; ++i) {
        os << boost::format("memuse: %4d %-55s %12d %12d %12d %12d\n")
          % i % vm_labels_[i]
          % vm_tmp_[i] % vm_ave_[i] % vm_min_[i] % vm_max_[i];
      }
#endif
      for (int i = 0; i < size; ++i) {
        if (counts_max_[i] > 0) {
          os << boost::format("timer: %5d %-55s %12.3lf %12.3lf %12.3lf %12.3lf"
                              " %10ld %10ld %10ld %10ld\n")
            % i % base_type::get_labels()[i]
            % base_type::get_measurements()[i] % sums_ave_[i] % sums_min_[i] % sums_max_[i]
            % base_type::get_counts()[i] % counts_ave_[i] % counts_min_[i] % counts_max_[i];
        }
      }
    }
  }
  void detailed_report(std::size_t interval = 1, std::ostream& os = std::clog) {
    #ifdef ALPS_ENABLE_TIMER_DETAILED
    int np, rank;
    MPI_Comm_size(comm_, &np);
    MPI_Comm_rank(comm_, &rank);
    if (rank == 0 && base_type::d_count_ == 0) {
      os << "timer: interval = " << interval << std::endl;
    }
    int size = base_type::d_counts_.size();
    if ((base_type::d_count_ + 1) % interval == 0) {
      if (rank == 0) {
        int buff_size = size * np;
        buff_counts_.resize(buff_size);
        buff_sums_.resize(buff_size);
      }
      MPI_Gather(&base_type::d_counts_[0], size, MPI_INT, &buff_counts_[0], size, MPI_INT, 0, comm_);
      MPI_Gather(&base_type::d_sums_[0], size, MPI_DOUBLE, &buff_sums_[0], size, MPI_DOUBLE, 0, comm_);
      if (rank == 0) {
        for (int i = 0; i < base_type::labels_.size(); ++i) {
          if (base_type::d_mapping_[i] >= 0) {
            for (int p = 0; p < np; ++p) {
              int k = size * p + base_type::d_mapping_[i];
              os << boost::format("detail: %d %d %d %.10f %d\n")
                % base_type::d_count_ % i % p % buff_sums_[k] % buff_counts_[k];
            }
          }
        }
      }
    }
    ++base_type::d_count_;
    std::fill(base_type::d_counts_.begin(), base_type::d_counts_.end(), 0);
    std::fill(base_type::d_sums_.begin(), base_type::d_sums_.end(), 0);
    MPI_Barrier(comm_);
    #endif
  }

private:
  MPI_Comm comm_;
  mutable std::valarray<double> counts_tmp_, counts_ave_, counts_min_, counts_max_,
    sums_tmp_, sums_ave_, sums_min_, sums_max_;
  mutable std::vector<int>  vm_tmp_, vm_ave_, vm_min_, vm_max_;
  #ifdef ALPS_ENABLE_TIMER_DETAILED
  mutable std::vector<int> buff_counts_;
  mutable std::vector<double> buff_sums_;
  #endif
  #ifdef ALPS_ENABLE_TIMER_BARRIER
  std::vector<int> barrier_;
  #endif
};

} // namespace detail

typedef detail::timer_mpi<detail::clock_mpi> timer_mpi;

} // namespace parapack
} // namespace alps

#else

namespace alps {
namespace parapack {

class timer_mpi : public timer {
public:
  BOOST_STATIC_CONSTANT(int, start_barrier = (1 << 1));
  BOOST_STATIC_CONSTANT(int, stop_barrier = (1 << 2));
  BOOST_STATIC_CONSTANT(int, barrier = (start_barrier | stop_barrier));
  timer_mpi(MPI_Comm, int = 1) {}
#ifndef ALPS_INDEP_SOURCE
  timer_mpi(MPI_Comm, alps::Parameters const&) {}
#endif
};

} // namespace parapack
} // namespace alps

#endif

#endif // ALPS_PARAPACK_TIMER_MPI_HPP
