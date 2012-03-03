/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>,
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

#ifndef LOOPER_PARALLEL_H
#define LOOPER_PARALLEL_H

#ifndef ALPS_INDEP_SOURCE
# include <alps/config.h>
# if defined(LOOPER_ENABLE_OPENMP) && defined(ALPS_ENABLE_OPENMP_WORKER) && !defined(LOOPER_OPENMP)
#  define LOOPER_OPENMP
# endif
#else
# if defined(LOOPER_ENABLE_OPENMP) && defined(_OPENMP) && !defined(LOOPER_OPENMP)
#  define LOOPER_OPENMP
# endif
#endif

// #define DEBUG_OUTPUT
// #define COMMUNICATION_TEST
// #define COMMUNICATION_DEBUG_OUTPUT

#include "expand.h"
#include "prime_factorization.h"
#include "union_find.h"
#include "timer_mpi.hpp"

#include <boost/detail/workaround.hpp>
#if !defined(BOOST_SPIRIT_USE_OLD_NAMESPACE)
# define BOOST_SPIRIT_USE_OLD_NAMESPACE
#endif
#include <boost/spirit/include/classic_actor.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_confix.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/throw_exception.hpp>
#include <numeric>
#include <stdexcept>
#include <mpi.h>

#ifdef LOOPER_OPENMP
# include <omp.h>
#endif

namespace looper {

namespace parallel {

template<typename LINK, typename COLLECTOR, typename ESTIMATE>
class chunk {
public:
  typedef std::pair<int, int> range_t;
  typedef LINK link_t;
  typedef COLLECTOR collector_t;
  typedef ESTIMATE estimate_t;

  chunk() {}
  explicit chunk(int ns) : links_(2 * ns), num_boundaries_(2 * ns) {
#ifndef COMMUNICATION_TEST
    coll_.clear_range();
#endif
  }
  void init(int ns) {
    links_.resize(2 * ns);
    num_boundaries_ = 2 * ns;
#ifndef COMMUNICATION_TEST
    coll_.clear_range();
#endif
  }
  static void init_timer(alps::parapack::timer_mpi& timer) {
#ifndef COMMUNICATION_TEST
    timer.registrate(51,  "      chunk::unify:0+1:all");
    timer.registrate(52,  "       chunk::unify:0+1:chunks_out.links_.resize");
    timer.registrate(53,  "       chunk::unify:0+1:chunks_out.prepare_links", timer.detailed);
    timer.registrate(54,  "       chunk::unify:0+1:unify_omp", timer.detailed);
    timer.registrate(55,  "       chunk::unify:0+1:set_root_omp", timer.detailed);
    timer.registrate(56,  "       chunk::unify:0+1:assign_cid", timer.detailed);
    timer.registrate(57,  "       chunk::unify:0+1:init-estimates_tg", timer.detailed);
    timer.registrate(58,  "       chunk::unify:0+1:collect_estimates", timer.detailed | timer.start_barrier);
    timer.registrate(59,  "       chunk::unify:0+1:unify_estimates", timer.detailed | timer.start_barrier);
    timer.registrate(60,  "       chunk::unify:0+1:construct_flip", timer.detailed | timer.start_barrier);
    timer.registrate(61,  "       chunk::unify:0+1:accumulate_measurements");
    timer.registrate(62,  "       chunk::unify:0+1:links_.resize");
    timer.registrate(63,  "       chunk::unify:0+1:estimates_.resize");
    timer.registrate(71,  "      chunk::unify:0+1+2:all");
    timer.registrate(81,  "      chunk::unify:0+1+2+3:all");
    timer.registrate(101, "      chunk::prepare_links");
    timer.registrate(111, "      chunk::assign_cid");
    timer.registrate(121, "      chunk::collect_estimates");
    timer.registrate(131, "      chunk::unify_estimates-all");
    timer.registrate(132, "       estimates_.resize(nc)");
    timer.registrate(133, "       for_estimates_");
    timer.registrate(141, "      chunk::construct_flip");
#endif // COMMUNICATION_TEST
  }

  int num_boundaries() const { return num_boundaries_; }

#ifndef COMMUNICATION_TEST
  template<typename FRAGMENT>
  void set(int pos, std::vector<FRAGMENT> const& fragments, int lower_offset,
           int upper_offset, collector_t const& coll, std::vector<estimate_t> const& estimates) {
    const int nb = num_boundaries();
    const int ns = nb / 2;
    const int noc = coll.num_open_clusters();
    if (lower_offset == 0 && upper_offset == ns) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < 2 * ns; ++v) {
        if (fragments[v].is_root()) {
          links_[v].set_id(fragments[v].id());
        } else {
          links_[v].set_parent(fragments[v].parent());
        }
      }
    } else if (lower_offset < upper_offset) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) {
        int vn = v;
        int vo = lower_offset + v;
        if (fragments[vo].is_root()) {
          links_[vn].set_id(fragments[vo].id());
        } else {
          int ro = root_index(fragments, vo);
          if (ro < upper_offset)
            links_[vn].set_parent(ro - lower_offset);
          else
            links_[vn].set_parent(ro + ns - upper_offset);
        }
      }
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) {
        int vn = ns + v;
        int vo = upper_offset + v;
        if (fragments[vo].is_root()) {
          links_[vn].set_id(fragments[vo].id());
        } else {
          int ro = root_index(fragments, vo);
          if (ro < upper_offset)
            links_[vn].set_parent(ro - lower_offset);
          else
            links_[vn].set_parent(ro + ns - upper_offset);
        }
      }
    } else { // (lower_offset > upper_offset)
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) {
        int vn = v;
        int vo = lower_offset + v;
        if (fragments[vo].is_root()) {
          links_[vn].set_id(fragments[vo].id());
        } else {
          int ro = root_index(fragments, vo);
          if (ro >= lower_offset)
            links_[vn].set_parent(ro - lower_offset);
          else
            links_[vn].set_parent(ro + ns - upper_offset);
        }
      }
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) {
        int vn = ns + v;
        int vo = upper_offset + v;
        if (fragments[vo].is_root()) {
          links_[vn].set_id(fragments[vo].id());
        } else {
          int ro = root_index(fragments, vo);
          if (ro >= lower_offset)
            links_[vn].set_parent(ro - lower_offset);
          else
            links_[vn].set_parent(ro + ns - upper_offset);
        }
      }
    }
    coll_ = coll;
    coll_.set_range(pos);
    looper::expand(estimates_, noc);
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int c = 0; c < noc; ++c) estimates_[c] = estimates[c];
  }

  bool empty() const { return coll_.empty(); }
  range_t const& range() const { return coll_.range(); }
  bool operator<(chunk const& rhs) const { return range().second + 1 < rhs.range().first; }
  bool operator<=(chunk const& rhs) const { return range().second + 1 == rhs.range().first; }

  void clear() { coll_.clear_range(); }
  void move_from(chunk& from) {
    links_.swap(from.links_);
    coll_ = from.coll_;
    estimates_.swap(from.estimates_);
  }
  void copy_from(chunk const& from) {
    const int nb = num_boundaries();
    const int noc = from.coll_.num_open_clusters();
    coll_ = from.coll_;
    looper::expand(estimates_, noc);
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int v = 0; v < nb; ++v) links_[v] = from.links_[v];
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int c = 0; c < noc; ++c) estimates_[c] = from.estimates_[c];
    }
  }
#endif // COMMUNICATION_TEST

  static void sendrecv(chunk& chunk_s, int dest, chunk& chunk_r, int source, int tag,
                       MPI_Comm comm) {
#ifdef COMMUNICATION_DEBUG_OUTPUT
    int pid;
    MPI_Comm_rank(comm, &pid);
    std::clog << pid << ": sendrecv (dest = " << dest << ", source = " << source << ")\n";
#endif
    MPI_Status status;
#ifndef COMMUNICATION_TEST
    MPI_Sendrecv(&chunk_s.coll_, sizeof(collector_t), MPI_BYTE, dest, tag,
                 &chunk_r.coll_, sizeof(collector_t), MPI_BYTE, source, tag,
                 comm, &status);
#endif // COMMUNICATION_TEST
    const int nb = chunk_s.num_boundaries();
    MPI_Sendrecv(&chunk_s.links_[0], nb, MPI_INT, dest, tag,
                 &chunk_r.links_[0], nb, MPI_INT, source, tag,
                 comm, &status);
#ifndef COMMUNICATION_TEST
    const int noc = chunk_r.coll_.num_open_clusters();
    looper::expand(chunk_r.estimates_, noc);
    const int count_s = sizeof(estimate_t) * chunk_s.coll_.num_open_clusters();
    const int count_r = sizeof(estimate_t) * chunk_r.coll_.num_open_clusters();
    MPI_Sendrecv(&(chunk_s.estimates_[0]), count_s, MPI_BYTE, dest, tag,
                 &(chunk_r.estimates_[0]), count_r, MPI_BYTE, source, tag,
                 comm, &status);
#endif // COMMUNICATION_TEST
  }
  static void sendrecv2(chunk& chunk_s0, int dest0, chunk& chunk_s1, int dest1,
                        chunk& chunk_r0, int source0, chunk& chunk_r1, int source1,
                        int tag, MPI_Comm comm) {
#ifdef COMMUNICATION_DEBUG_OUTPUT
    int pid;
    MPI_Comm_rank(comm, &pid);
    std::clog << pid << ": sendrecv2 (dest0 = " << dest0 << ", dest1 = " << dest1
              << ", source0 = " << source0 << ", source1 = " << source1 << ")\n";
#endif
    MPI_Status status[2];
#ifndef COMMUNICATION_TEST
    MPI_Request send0[2], recv0[2];
    MPI_Isend(&chunk_s0.coll_, sizeof(collector_t), MPI_BYTE, dest0, tag, comm, &send0[0]);
    MPI_Isend(&chunk_s1.coll_, sizeof(collector_t), MPI_BYTE, dest1, tag, comm, &send0[1]);
    MPI_Irecv(&chunk_r0.coll_, sizeof(collector_t), MPI_BYTE, source0, tag, comm, &recv0[0]);
    MPI_Irecv(&chunk_r1.coll_, sizeof(collector_t), MPI_BYTE, source1, tag, comm, &recv0[1]);
#endif // COMMUNICATION_TEST
    const int nb = chunk_s0.num_boundaries();
    MPI_Request send1[2], recv1[2];
    MPI_Isend(&chunk_s0.links_[0], nb, MPI_INT, dest0, tag, comm, &send1[0]);
    MPI_Isend(&chunk_s1.links_[0], nb, MPI_INT, dest1, tag, comm, &send1[1]);
    MPI_Irecv(&chunk_r0.links_[0], nb, MPI_INT, source0, tag, comm, &recv1[0]);
    MPI_Irecv(&chunk_r1.links_[0], nb, MPI_INT, source1, tag, comm, &recv1[1]);
#ifndef COMMUNICATION_TEST
    MPI_Waitall(2, send0, status);
    MPI_Waitall(2, recv0, status);
    looper::expand(chunk_r0.estimates_, chunk_r0.coll_.num_open_clusters());
    looper::expand(chunk_r1.estimates_, chunk_r1.coll_.num_open_clusters());
    const int count_s0 = sizeof(estimate_t) * chunk_s0.coll_.num_open_clusters();
    const int count_s1 = sizeof(estimate_t) * chunk_s1.coll_.num_open_clusters();
    const int count_r0 = sizeof(estimate_t) * chunk_r0.coll_.num_open_clusters();
    const int count_r1 = sizeof(estimate_t) * chunk_r1.coll_.num_open_clusters();
    MPI_Request send2[2], recv2[2];
    MPI_Isend(&(chunk_s0.estimates_[0]), count_s0, MPI_BYTE, dest0, tag, comm, &send2[0]);
    MPI_Isend(&(chunk_s1.estimates_[0]), count_s1, MPI_BYTE, dest1, tag, comm, &send2[1]);
    MPI_Irecv(&(chunk_r0.estimates_[0]), count_r0, MPI_BYTE, source0, tag, comm, &recv2[0]);
    MPI_Irecv(&(chunk_r1.estimates_[0]), count_r1, MPI_BYTE, source1, tag, comm, &recv2[1]);
#endif // COMMUNICATION_TEST
    MPI_Waitall(2, send1, status);
#ifndef COMMUNICATION_TEST
    MPI_Waitall(2, send2, status);
#endif // COMMUNICATION_TEST
    MPI_Waitall(2, recv1, status);
#ifndef COMMUNICATION_TEST
    MPI_Waitall(2, recv2, status);
#endif // COMMUNICATION_TEST
  }

#ifndef COMMUNICATION_TEST
  collector_t const& get_collector() const { return coll_; }

  void prepare_links(const std::vector<link_t>& links_in, int start_l, int start_u,
    alps::parapack::timer_mpi& timer) {
    timer.start(101);
    const int nb = num_boundaries();
    const int ns = nb / 2;
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) {
        links_[v + start_l] = links_in[v];
        if (!links_in[v].is_root()) {
          if (links_in[v].parent() >= ns)
            links_[v + start_l].set_parent(links_in[v].parent() + start_u - ns);
          else
            links_[v + start_l].set_parent(links_in[v].parent() + start_l);
        }
      }
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int v = ns; v < 2 * ns; ++v) {
        links_[v + start_u - ns] = links_in[v];
        if (!links_in[v].is_root()) {
          if (links_in[v].parent() >= ns)
            links_[v + start_u - ns].set_parent(links_in[v].parent() + start_u - ns);
          else
            links_[v + start_u - ns].set_parent(links_in[v].parent() + start_l);
        }
      }
    }
    timer.stop(101);
  }

  int assign_cid(int start, int end, int nc_in, alps::parapack::timer_mpi& timer) {
    timer.start(111);
    #ifdef LOOPER_OPENMP
    int num_threads = omp_get_max_threads();
    std::vector<int> ncl(num_threads);
    for (int tid = 0; tid < num_threads; ++tid) ncl[tid] = 0;

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int total = 0;
      int offset;
      #pragma omp for schedule(static) nowait
      for (int v = start; v < end; ++v)
        total += links_[v].is_root();

      ncl[tid] = total;
      #pragma omp barrier
      offset = nc_in;
      for (int i = 0; i < tid; ++i) offset += ncl[i];

      #pragma omp for schedule(static)
      for (int v = start; v < end; ++v) {
        int mask = links_[v].is_root();
        links_[v].set_id((mask * (offset)) + (1 - mask) * links_[v].id());
        offset += mask;
      }
      if (tid == num_threads - 1) nc_in = offset;
    }
#else
    for (int v = start; v < end; ++v) if (links_[v].is_root()) links_[v].set_id(nc_in++);
#endif
    timer.stop(111);
    return nc_in;
  }

  // should be run in parallel region
  void collect_estimates(chunk const& chunk_in, int start, int end, int offset,
    std::vector<estimate_t>& estimates) {
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int v = start; v < end; ++v) {
      if (chunk_in.links_[v - offset].is_root()) {
        int cid_in = chunk_in.links_[v - offset].id();
        int cid_out = cluster_id(links_, v);
        estimates[cid_out] += chunk_in.estimates_[cid_in];
      }
    }
  }

#ifdef LOOPER_OPENMP
  void unify_estimates(int nc, std::vector<std::vector<estimate_t> > const& estimates_tg,
    alps::parapack::timer_mpi& timer) {
    timer.start(131);
    const int num_threads = omp_get_max_threads();
    timer.start(132);
    looper::expand(estimates_, nc);
    timer.stop(132);
    timer.start(133);
    #pragma omp parallel for schedule(static)
    for (int c = 0; c < nc; ++c)
      estimates_[c] = estimates_tg[0][c];
    for (int p = 1; p < num_threads; ++p) {
      #pragma omp parallel for schedule(static)
      for (int c = 0; c < nc; ++c)
        estimates_[c] += estimates_tg[p][c];
    }
    timer.stop(133);
    timer.stop(131);
  }
#endif

  template<typename FLIP>
  void construct_flip(const std::vector<link_t>& links_in, std::vector<FLIP>& flip_table,
    std::vector<int>& map_table, int noc, int start, int end, int offset,
    alps::parapack::timer_mpi& timer) {
    timer.start(141);
    if (map_table.size() > 0) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = start; v < end; ++v) {
        int vo = v - offset;
        if (links_in[vo].is_root()) {
          int new_id = cluster_id(links_, v);
          int old_id = links_in[vo].id();
          if (new_id < noc)
            flip_table[old_id].set_local_cid(new_id);
          else
            flip_table[old_id].set_flip(estimates_[new_id].to_flip);
          map_table[old_id] = new_id;
        }
      }
    } else {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = start; v < end; ++v) {
        int vo = v - offset;
        if (links_in[vo].is_root()) {
          int new_id = cluster_id(links_, v);
          int old_id = links_in[vo].id();
          if (new_id < noc)
            flip_table[old_id].set_local_cid(new_id);
          else
            flip_table[old_id].set_flip(estimates_[new_id].to_flip);
        }
      }
    }
    timer.stop(141);
  }

  #ifdef LOOPER_OPENMP
  template<typename FLIP>
  static void unify(chunk const& chunk_in0, chunk const& chunk_in1, bool connect_periodic,
    int flip_index, chunk& chunk_out, std::vector<FLIP>& flip_table, std::vector<int>& map_table, 
    alps::parapack::timer_mpi& timer, std::vector<std::vector<estimate_t> >& estimates_tg) {
  #else
  template<typename FLIP>
  static void unify(chunk const& chunk_in0, chunk const& chunk_in1, bool connect_periodic,
    int flip_index, chunk& chunk_out, std::vector<FLIP>& flip_table, std::vector<int>& map_table, 
    alps::parapack::timer_mpi& timer) {
  #endif
    timer.start(51);
    const int nb = chunk_in0.num_boundaries();
    const int ns = nb / 2; // (nb should be even)

    // prepare chunk_out.links_
    // [0L|0U] [1L|1U] -> [0L|1U|1L|0U]
    timer.start(52);
    looper::expand(chunk_out.links_, 2 * nb);
    timer.stop(52);
    timer.start(53);
    chunk_out.prepare_links(chunk_in0.links_, 0, 3 * ns, timer);      // for [0L|0U]
    chunk_out.prepare_links(chunk_in1.links_, 2 * ns, 1 * ns, timer); // for [1L|1U]
    timer.stop(53);

    // connect between lower and upper part
    timer.start(54);
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int v = 0; v < ns; ++v)
      looper::union_find::unify(chunk_out.links_, 2 * ns + v, 3 * ns + v);
    if (connect_periodic) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) looper::union_find::unify(chunk_out.links_, v, ns + v);
    }
    timer.stop(54);
    timer.start(55);
    pack_tree(chunk_out.links_, 2 * ns);
    timer.stop(55);

    // assign cluster ID
    timer.start(56);
    int nc = 0; // total number of clusters
    nc = chunk_out.assign_cid(0, 2 * ns, 0, timer);
    int noc = connect_periodic ? 0 : nc; // number of open clusters
    nc = chunk_out.assign_cid(2 * ns, 4 * ns, nc, timer);
    timer.stop(56);

    // unify estimates
    timer.start(57);
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      std::vector<estimate_t>& estimates = estimates_tg[omp_get_thread_num()];
      #else
      std::vector<estimate_t>& estimates = chunk_out.estimates_;
      #endif
      looper::expand(estimates, nc);
      estimate_t estimate_init;
      for (int c = 0; c < nc; ++c) estimates[c] = estimate_init;
    }
    timer.stop(57);
    timer.start(58);
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      timer.prof_start();
      #ifdef LOOPER_OPENMP
      std::vector<estimate_t>& estimates = estimates_tg[omp_get_thread_num()];
      #else
      std::vector<estimate_t>& estimates = chunk_out.estimates_;
      #endif
      chunk_out.collect_estimates(chunk_in0, 0, ns, 0, estimates); // 0L
      chunk_out.collect_estimates(chunk_in1, ns, 2 * ns, 0, estimates); // 1U
      chunk_out.collect_estimates(chunk_in1, 2 * ns, 3 * ns, 2 * ns, estimates); // 1L
      chunk_out.collect_estimates(chunk_in0, 3 * ns, 4 * ns, 2 * ns, estimates); // 0U
      timer.prof_stop();
    }
    timer.stop(58);
    timer.start(59);
    #ifdef LOOPER_OPENMP
    timer.prof_start();
    chunk_out.unify_estimates(nc, estimates_tg, timer);
    timer.prof_stop();
    #endif
    timer.stop(59);

    // construct flip table
    timer.start(60);
    timer.prof_start();
    if (flip_index == 0) {
      // for 0L and 0U
      chunk_out.construct_flip(chunk_in0.links_, flip_table, map_table, noc, 0, ns, 0, timer);
      chunk_out.construct_flip(chunk_in0.links_, flip_table, map_table, noc, 3 * ns, 4 * ns,
                               2 * ns, timer);
    } else if (flip_index == 1) {
      // for 1L and 1U
      chunk_out.construct_flip(chunk_in1.links_, flip_table, map_table, noc, 2 * ns, 3 * ns,
                               2 * ns, timer);
      chunk_out.construct_flip(chunk_in1.links_, flip_table, map_table, noc, ns, 2 * ns, 0, timer);
    } else if (flip_index != -1) {
      std::cerr << "invalid flip_index " << flip_index << " in unify\n";
      boost::throw_exception(std::logic_error("chunk::unify"));
    }
    timer.prof_stop();
    timer.stop(60);

    // accumulate measurements of closed clusters
    timer.start(61);
    chunk_out.coll_ = chunk_in0.coll_;
    chunk_out.coll_ += chunk_in1.coll_;
#ifdef LOOPER_OPENMP
    int num_threads = omp_get_max_threads();
    std::vector<collector_t> coll_tmp(num_threads);
    #pragma omp parallel
    {
      collector_t coll;
      coll.set_range(chunk_out.coll_);
      #pragma omp for schedule(static)
      for (int c = noc; c < nc; ++c) coll += chunk_out.estimates_[c];
      coll_tmp[omp_get_thread_num()] = coll;
    }
    for (int tid = 0; tid < num_threads; ++tid) chunk_out.coll_ += coll_tmp[tid];
#else
    for (int c = noc; c < nc; ++c) chunk_out.coll_ += chunk_out.estimates_[c];
#endif
    chunk_out.coll_.set_num_open_clusters(noc);
    chunk_out.coll_.inc_num_clusters(nc - noc);
    timer.stop(61);
    timer.start(62);
    timer.stop(62);
    timer.start(63);
    // if (chunk_out.estimates_.size() < noc) chunk_out.estimates_.resize(noc); // not necessary?
    timer.stop(63);
    timer.stop(51);
  }

#ifdef LOOPER_OPENMP
  template<typename FLIP>
  static void unify(chunk const& chunk_in0, chunk const& chunk_in1, chunk const& chunk_in2,
    bool connect_periodic, int flip_index, chunk& chunk_out,  std::vector<FLIP>& flip_table,
    std::vector<int>& map_table,                   
    alps::parapack::timer_mpi& timer, std::vector<std::vector<estimate_t> >& estimates_tg) {
#else
  template<typename FLIP>
  static void unify(chunk const& chunk_in0, chunk const& chunk_in1, chunk const& chunk_in2,
    bool connect_periodic, int flip_index, chunk& chunk_out,  std::vector<FLIP>& flip_table,
    std::vector<int>& map_table,                   
    alps::parapack::timer_mpi& timer) {
#endif
    timer.start(71);
    const int nb = chunk_in0.num_boundaries();
    const int ns = nb / 2; // (nb should be even)

    // prepare chunk_out.links_
    // [0L|0U] [1L|1U] [2L|2U] -> [0L|2U|2L|0U|1L|1U]
    looper::expand(chunk_out.links_, 3 * nb);
    chunk_out.prepare_links(chunk_in0.links_, 0, 3 * ns, timer);      // for [0L|0U]
    chunk_out.prepare_links(chunk_in1.links_, 4 * ns, 5 * ns, timer); // for [1L|1U]
    chunk_out.prepare_links(chunk_in2.links_, 2 * ns, ns, timer);     // for [2L|2U]

    // connect between lower and upper part
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int v = 0; v < ns; ++v) {
      looper::union_find::unify(chunk_out.links_, 2 * ns + v, 5 * ns + v);
      looper::union_find::unify(chunk_out.links_, 3 * ns + v, 4 * ns + v);
    }
    if (connect_periodic) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) looper::union_find::unify(chunk_out.links_, v, ns + v);
    }
    pack_tree(chunk_out.links_, 2 * ns);

    // assign cluster ID
    int nc = 0; // total number of clusters
    nc = chunk_out.assign_cid(0, 2 * ns, 0, timer);
    int noc = connect_periodic ? 0 : nc; // number of open clusters
    nc = chunk_out.assign_cid(2 * ns, 6 * ns, nc, timer);

    // unify estimates
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      std::vector<estimate_t>& estimates = estimates_tg[omp_get_thread_num()];
      #else
      std::vector<estimate_t>& estimates = chunk_out.estimates_;
      #endif
      looper::expand(estimates, nc);
      estimate_t estimate_init;
      for (int c = 0; c < nc; ++c) estimates[c] = estimate_init;
    }
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      std::vector<estimate_t>& estimates = estimates_tg[omp_get_thread_num()];
      #else
      std::vector<estimate_t>& estimates = chunk_out.estimates_;
      #endif
      chunk_out.collect_estimates(chunk_in0, 0, ns, 0, estimates); // 0L
      chunk_out.collect_estimates(chunk_in2, ns, 2 * ns, 0, estimates); // 2U
      chunk_out.collect_estimates(chunk_in2, 2 * ns, 3 * ns, 2 * ns, estimates); // 2L
      chunk_out.collect_estimates(chunk_in0, 3 * ns, 4 * ns, 2 * ns, estimates); // 0U
      chunk_out.collect_estimates(chunk_in1, 4 * ns, 6 * ns, 4 * ns, estimates); // 1L & 1U
    }
    #ifdef LOOPER_OPENMP
    chunk_out.unify_estimates(nc, estimates_tg, timer);
    #endif

    // construct flip table
    if (flip_index == 0) {
      // 0L and 0U
      chunk_out.construct_flip(chunk_in0.links_, flip_table, map_table, noc, 0, ns, 0, timer);
      chunk_out.construct_flip(chunk_in0.links_, flip_table, map_table, noc, 3 * ns, 4 * ns,
                               2 * ns, timer);
    } else if (flip_index == 1) {
      // 1L and 1U
      chunk_out.construct_flip(chunk_in1.links_, flip_table, map_table, noc, 4 * ns, 6 * ns,
                               4 * ns, timer);
    } else if (flip_index == 2) {
      // 2L and 2U
      chunk_out.construct_flip(chunk_in2.links_, flip_table, map_table, noc, 2 * ns, 3 * ns,
                               2 * ns, timer);
      chunk_out.construct_flip(chunk_in2.links_, flip_table, map_table, noc, ns, 2 * ns, 0, timer);
    } else if (flip_index != -1) {
      std::cerr << "invalid flip_index " << flip_index << " in unify\n";
      boost::throw_exception(std::logic_error("chunk::unify"));
    }

    // accumulate measurements of closed clusters
    chunk_out.coll_ = chunk_in0.coll_;
    chunk_out.coll_ += chunk_in1.coll_;
    chunk_out.coll_ += chunk_in2.coll_;
#ifdef LOOPER_OPENMP
    int num_threads = omp_get_max_threads();
    std::vector<collector_t> coll_tmp(num_threads);
    #pragma omp parallel
    {
      collector_t coll;
      coll.set_range(chunk_out.coll_);
      #pragma omp for schedule(static)
      for (int c = noc; c < nc; ++c) coll += chunk_out.estimates_[c];
      coll_tmp[omp_get_thread_num()] = coll;
    }
    for (int tid = 0; tid < num_threads; ++tid) chunk_out.coll_ += coll_tmp[tid];
#else
    for (int c = noc; c < nc; ++c) chunk_out.coll_ += chunk_out.estimates_[c];
#endif
    chunk_out.coll_.set_num_open_clusters(noc);
    chunk_out.coll_.inc_num_clusters(nc - noc);
    // if (chunk_out.estimates_.size() < noc) chunk_out.estimates_.resize(noc); // not necessary?
    timer.stop(71);
  }

  #ifdef LOOPER_OPENMP
  template<typename FLIP>
  static void unify(chunk const& chunk_in0, chunk const& chunk_in1, chunk const& chunk_in2,
    chunk const& chunk_in3, bool connect_periodic, int flip_index, chunk& chunk_out,
    std::vector<FLIP>& flip_table, std::vector<int>& map_table, alps::parapack::timer_mpi& timer,
    std::vector<std::vector<estimate_t> >& estimates_tg) {
  #else
  template<typename FLIP>
  static void unify(chunk const& chunk_in0, chunk const& chunk_in1, chunk const& chunk_in2,
    chunk const& chunk_in3, bool connect_periodic, int flip_index, chunk& chunk_out,
    std::vector<FLIP>& flip_table, std::vector<int>& map_table, alps::parapack::timer_mpi& timer) {
  #endif
    timer.start(81);
    const int nb = chunk_in0.num_boundaries();
    const int ns = nb / 2; // (nb should be even)

    // prepare chunk_out.links_
    // [0L|0U] [1L|1U] [2L|2U] [3L|3U] -> [0L|3U|3L|0U|1L|1U|2L|2U]
    looper::expand(chunk_out.links_, 4 * nb);
    chunk_out.prepare_links(chunk_in0.links_, 0, 3 * ns, timer);      // for [0L|0U]
    chunk_out.prepare_links(chunk_in1.links_, 4 * ns, 5 * ns, timer); // for [1L|1U]
    chunk_out.prepare_links(chunk_in2.links_, 6 * ns, 7 * ns, timer); // for [2L|2U]
    chunk_out.prepare_links(chunk_in3.links_, 2 * ns, 1 * ns, timer); // for [3L|3U]

    // connect between lower and upper part
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int v = 0; v < ns; ++v) {
      looper::union_find::unify(chunk_out.links_, 2 * ns + v, 7 * ns + v);
      looper::union_find::unify(chunk_out.links_, 3 * ns + v, 4 * ns + v);
      looper::union_find::unify(chunk_out.links_, 5 * ns + v, 6 * ns + v);
    }
    if (connect_periodic) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int v = 0; v < ns; ++v) looper::union_find::unify(chunk_out.links_, v, ns + v);
    }
    pack_tree(chunk_out.links_, 2 * ns);

    // assign cluster ID
    int nc = 0; // total number of clusters
    nc = chunk_out.assign_cid(0, 2 * ns, 0, timer);
    int noc = connect_periodic ? 0 : nc; // number of open clusters
    nc = chunk_out.assign_cid(2 * ns, 8 * ns, nc, timer);

    // unify estimates
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      std::vector<estimate_t>& estimates = estimates_tg[omp_get_thread_num()];
      #else
      std::vector<estimate_t>& estimates = chunk_out.estimates_;
      #endif
      looper::expand(estimates, nc);
      estimate_t estimate_init;
      for (int c = 0; c < nc; ++c) estimates[c] = estimate_init;
    }
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      std::vector<estimate_t>& estimates = estimates_tg[omp_get_thread_num()];
      #else
      std::vector<estimate_t>& estimates = chunk_out.estimates_;
      #endif
      chunk_out.collect_estimates(chunk_in0, 0, ns, 0, estimates); // 0L
      chunk_out.collect_estimates(chunk_in3, ns, 2 * ns, 0, estimates); // 3L
      chunk_out.collect_estimates(chunk_in3, 2 * ns, 3 * ns, 2 * ns, estimates); // 3L
      chunk_out.collect_estimates(chunk_in0, 3 * ns, 4 * ns, 2 * ns, estimates); // 0U
      chunk_out.collect_estimates(chunk_in1, 4 * ns, 6 * ns, 4 * ns, estimates); // 1L & 1U
      chunk_out.collect_estimates(chunk_in2, 6 * ns, 8 * ns, 6 * ns, estimates); // 2L & 2U
    }
    #ifdef LOOPER_OPENMP
    chunk_out.unify_estimates(nc, estimates_tg, timer);
    #endif

    // construct flip table
    if (flip_index == 0) {
      // 0L and 0U
      chunk_out.construct_flip(chunk_in0.links_, flip_table, map_table, noc, 0, ns, 0, timer);
      chunk_out.construct_flip(chunk_in0.links_, flip_table, map_table, noc, 3 * ns, 4 * ns,
                               2 * ns, timer);
    } else if (flip_index == 3) {
      // 3L and 3U
      chunk_out.construct_flip(chunk_in3.links_, flip_table, map_table, noc, 2 * ns, 3 * ns,
                               2 * ns, timer);
      chunk_out.construct_flip(chunk_in3.links_, flip_table, map_table, noc, ns, 2 * ns, 0, timer);
    } else if (flip_index != -1) {
      std::cerr << "invalid flip_index " << flip_index << " in unify\n";
      boost::throw_exception(std::logic_error("chunk::unify"));
    }

    // accumulate measurements of closed clusters
    chunk_out.coll_ = chunk_in0.coll_;
    chunk_out.coll_ += chunk_in1.coll_;
    chunk_out.coll_ += chunk_in2.coll_;
    chunk_out.coll_ += chunk_in3.coll_;
#ifdef LOOPER_OPENMP
    int num_threads = omp_get_max_threads();
    std::vector<collector_t> coll_tmp(num_threads);
    #pragma omp parallel
    {
      collector_t coll;
      coll.set_range(chunk_out.coll_);
      #pragma omp for schedule(static)
      for (int c = noc; c < nc; ++c) coll += chunk_out.estimates_[c];
      coll_tmp[omp_get_thread_num()] = coll;
    }
    for (int tid = 0; tid < num_threads; ++tid) chunk_out.coll_ += coll_tmp[tid];
#else
    for (int c = noc; c < nc; ++c) chunk_out.coll_ += chunk_out.estimates_[c];
#endif
    chunk_out.coll_.set_num_open_clusters(noc);
    chunk_out.coll_.inc_num_clusters(nc - noc);
    // if (chunk_out.estimates_.size() < noc) chunk_out.estimates_.resize(noc); // not necessary?
    timer.stop(81);
  }

  // generate identity flip table
  template<typename FLIP>
    void fill_flip(std::vector<FLIP>& flip_table, std::vector<int>& map_table) const {
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int c = 0; c < coll_.num_open_clusters(); ++c) flip_table[c].set_local_cid(c);
    if (map_table.size() > 0) {
      #ifdef LOOPER_OPENMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int c = 0; c < coll_.num_open_clusters(); ++c) map_table[c] = c;
    }
  }

  void output_links(int pid, int stage, int step, std::string title, std::ostream& os) const {
    if (stage >= 0) {
      if (step >= 0) {
        title = std::string("stage ") + boost::lexical_cast<std::string>(stage) + ", step " +
          boost::lexical_cast<std::string>(step) + ": " + title;
      } else {
        title = std::string("stage ") + boost::lexical_cast<std::string>(stage) + ": " + title;
      }
    }
    os << "process " << pid << ": " << title << " [";
    for (int i = 0; i < num_boundaries(); ++i) {
      if (i != 0) os << ',';
      if (links_[i].is_root())
        os << '*';
      else
        os << links_[i].parent();
    }
    os << "] {";
    for (int i = 0; i < num_boundaries(); ++i) {
      if (i != 0) os << ',';
      if (links_[i].is_root())
        os << '(' << cluster_id(links_, i) << ')';
      else
        os << ' ' << cluster_id(links_, i) << ' ';
    }
    os << "}" << std::endl;
  }

  void output_collector(int pid, int stage, int step, std::string title, std::ostream& os) const {
    if (stage >= 0) {
      if (step >= 0) {
        title = std::string("stage ") + boost::lexical_cast<std::string>(stage) + ", step " +
          boost::lexical_cast<std::string>(step) + ": " + title;
      } else {
        title = std::string("stage ") + boost::lexical_cast<std::string>(stage) + ": " + title;
      }
    }
    os << "process " << pid << ": " << title << " range = [" << coll_.range().first << ":"
       << coll_.range().second << "], number of closed clusters = "
       << coll_.num_clusters() << ", number of open clusters = " << coll_.num_open_clusters()
       << std::endl;
  }
#endif // COMMUNICATION_TEST

private:
  std::vector<link_t> links_;
#ifndef COMMUNICATION_TEST
  collector_t coll_;
  std::vector<estimate_t> estimates_;
#endif // COMMUNICATION_TEST
  int num_boundaries_; // number of sites at imaginary-boundaries (2 * num_sites_)
};


template<typename LINK, typename COLLECTOR, typename ESTIMATE>
class chunks {
public:
  typedef LINK link_t;
  typedef COLLECTOR collector_t;
  typedef ESTIMATE estimate_t;
  typedef chunk<link_t, collector_t, estimate_t> chunk_t;

  chunks() {}
  explicit chunks(int ns) : chunk0_(ns), chunk1_(ns), chunk_tmp_(ns) {}
  void init(int ns) {
    chunk0_.init(ns);
    chunk1_.init(ns);
    chunk_tmp_.init(ns);
  }
  static void init_timer(alps::parapack::timer_mpi& timer) {
#ifndef COMMUNICATION_TEST
    timer.registrate(41, "     chunks::insert:all");
    timer.registrate(42, "      chunks::insert:1+1");
    timer.registrate(43, "      chunks::insert:1+2");
    timer.registrate(44, "      chunks::insert:2+1");
    timer.registrate(45, "      chunks::insert:2+2");
#endif // COMMUNICATION_TEST
    chunk_t::init_timer(timer);
  }
#ifndef COMMUNICATION_TEST
  template<typename FRAGMENT>
  void set(int pos, std::vector<FRAGMENT> const& fragments, int lower_offset,
    int upper_offset, collector_t const& coll, std::vector<estimate_t> const& estimates) {
    orig_ = pos;
    chunk0_.set(pos, fragments, lower_offset, upper_offset, coll, estimates);
    chunk1_.clear();
  }
  bool empty() const { return chunk0_.empty(); }
  int size() const { return chunk0_.empty() ? 0 : (chunk1_.empty() ? 1 : 2); }
  void clear() { chunk0_.clear(); chunk1_.clear(); }
  void move_from(chunks& from) {
    orig_ = from.orig_;
    chunk0_.move_from(from.chunk0_);
    from.chunk0_.clear();
    if (!from.chunk1_.empty()) {
      chunk1_.move_from(from.chunk1_);
      from.chunk1_.clear();
    } else {
      chunk1_.clear();
    }
  }
  void copy_from(chunks const& from) {
    orig_ = from.orig_;
    chunk0_.copy_from(from.chunk0_);
    if (!from.chunk1_.empty()) {
      chunk1_.copy_from(from.chunk1_);
    } else {
      chunk1_.clear();
    }
  }
#endif // COMMUNICATION_TEST

  static void sendrecv(chunks& chunks_s, int dest, chunks& chunks_r, int source, int tag,
                       MPI_Comm comm) {
    chunk_t::sendrecv(chunks_s.chunk0_, dest, chunks_r.chunk0_, source, tag, comm);
#ifndef COMMUNICATION_TEST
    chunks_r.chunk1_.clear();
#endif // COMMUNICATION_TEST
  }
  static void sendrecv2(chunks& chunks_s, int dest, chunks& chunks_r, int source, int tag,
                        MPI_Comm comm) {
#ifndef COMMUNICATION_TEST
    if (chunks_s.size() == 1) {
      chunk_t::sendrecv2(chunks_s.chunk0_, dest, chunks_s.chunk0_, source,
                         chunks_r.chunk0_, source, chunks_r.chunk1_, dest, tag, comm);
    } else {
      chunk_t::sendrecv2(chunks_s.chunk0_, dest, chunks_s.chunk1_, source,
                         chunks_r.chunk0_, source, chunks_r.chunk1_, dest, tag, comm);
    }
#else
    chunk_t::sendrecv2(chunks_s.chunk0_, dest, chunks_s.chunk1_, source,
                       chunks_r.chunk0_, source, chunks_r.chunk1_, dest, tag, comm);
#endif // COMMUNICATION_TEST
  }

#ifndef COMMUNICATION_TEST
  collector_t const& get_collector() const {
    if (size() != 1) std::cerr << "get_collector: logic error " << size() << std::endl;
    return chunk0_.get_collector();
  }

  int get_cluster_id(int c) const { return chunk0_.get_cluster_id(c); }

  #ifdef LOOPER_OPENMP
  template<typename FLIP>
  void insert(chunks const& chunks_in, bool connect_periodic, std::vector<FLIP>& flip_table,
    std::vector<int>& map_table,
    alps::parapack::timer_mpi& timer, std::vector<std::vector<estimate_t> >& estimates_tg) {
  #else
  template<typename FLIP>
  void insert(chunks const& chunks_in, bool connect_periodic, std::vector<FLIP>& flip_table,
    std::vector<int>& map_table,
    alps::parapack::timer_mpi& timer) {
  #endif
    timer.start(41);
    bool valid = true;
    if (this->size() == 1) {
      if (chunks_in.size() == 1) {
        // 1 + 1
        timer.start(42);
        if (chunk0_ <= chunks_in.chunk0_) {
          // [old0][new0]
          #ifdef LOOPER_OPENMP
          chunk_t::unify(chunk0_, chunks_in.chunk0_, connect_periodic, 0, chunk_tmp_, flip_table,
                         map_table, timer, estimates_tg);
          #else
          chunk_t::unify(chunk0_, chunks_in.chunk0_, connect_periodic, 0, chunk_tmp_, flip_table,
                         map_table, timer);
          #endif
          chunk0_.move_from(chunk_tmp_);
        } else if (chunk0_ < chunks_in.chunk0_) {
          // [old0] ... [new0]
          chunk1_.copy_from(chunks_in.chunk0_);
          chunk0_.fill_flip(flip_table, map_table);
        } else if (chunks_in.chunk0_ <= chunk0_) {
          // [new0][old0]
          #ifdef LOOPER_OPENMP
          chunk_t::unify(chunks_in.chunk0_, chunk0_, connect_periodic, 1, chunk_tmp_, flip_table,
                         map_table, timer, estimates_tg);
          #else
          chunk_t::unify(chunks_in.chunk0_, chunk0_, connect_periodic, 1, chunk_tmp_, flip_table,
                         map_table, timer);
          #endif
          chunk0_.move_from(chunk_tmp_);
        } else if (chunks_in.chunk0_ < chunk0_) {
          // [new0] ... [old0]
          chunk1_.move_from(chunk0_);
          chunk0_.copy_from(chunks_in.chunk0_);
          chunk1_.fill_flip(flip_table, map_table);
        } else {
          valid = false;
        }
        timer.stop(42);
      } else {
        // 1 + 2
        timer.start(43);
        if (chunk0_ <= chunks_in.chunk0_) {
          if (chunks_in.chunk0_ <= chunks_in.chunk1_) {
            // [old0][new0][new1]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunk0_, chunks_in.chunk0_, chunks_in.chunk1_, connect_periodic, 0,
                           chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunk0_, chunks_in.chunk0_, chunks_in.chunk1_, connect_periodic, 0,
                           chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
          } else if (chunks_in.chunk0_ < chunks_in.chunk1_) {
            // [old0][new0] ... [new1]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunk0_, chunks_in.chunk0_, 0, 0, chunk_tmp_, flip_table, map_table,
                           timer, estimates_tg);
            #else
            chunk_t::unify(chunk0_, chunks_in.chunk0_, 0, 0, chunk_tmp_, flip_table, map_table,
                           timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
            chunk1_.copy_from(chunks_in.chunk1_);
          } else if (chunks_in.chunk1_ <= chunk0_) {
            // [new1][old0][new0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk1_, chunk0_, chunks_in.chunk0_, connect_periodic, 1,
                           chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk1_, chunk0_, chunks_in.chunk0_, connect_periodic, 1,
                           chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
          } else {
            valid = false;
          }
        } else if (chunk0_ <= chunks_in.chunk1_) {
          if (chunks_in.chunk1_ <= chunks_in.chunk0_) {
            // [old0][new1][new0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunk0_, chunks_in.chunk1_, chunks_in.chunk0_, connect_periodic, 0,
                           chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunk0_, chunks_in.chunk1_, chunks_in.chunk0_, connect_periodic, 0,
                           chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
          } else if (chunks_in.chunk1_ < chunks_in.chunk0_) {
            // [old0][new1] ... [new0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunk0_, chunks_in.chunk1_, 0, 0, chunk_tmp_, flip_table, map_table,
                           timer, estimates_tg);
            #else
            chunk_t::unify(chunk0_, chunks_in.chunk1_, 0, 0, chunk_tmp_, flip_table, map_table,
                           timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
            chunk1_.copy_from(chunks_in.chunk0_);
          } else if (chunks_in.chunk1_ <= chunk0_) {
            // [new0][old0][new1]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk0_, chunk0_, chunks_in.chunk1_,connect_periodic, 1,
                           chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk0_, chunk0_, chunks_in.chunk1_,connect_periodic, 1,
                           chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
          } else {
            valid = false;
          }
        } else if (chunks_in.chunk0_ <= chunk0_) {
          if (chunks_in.chunk1_ <= chunks_in.chunk0_) {
            // [new1][new0][old0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk1_, chunks_in.chunk0_, chunk0_, connect_periodic, 2,
                           chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk1_, chunks_in.chunk0_, chunk0_, connect_periodic, 2,
                           chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
          } else if (chunks_in.chunk1_ < chunks_in.chunk0_) {
            // [new1] ... [new0][old0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk0_, chunk0_, connect_periodic, 1, chunk_tmp_, flip_table,
                           map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk0_, chunk0_, connect_periodic, 1, chunk_tmp_, flip_table,
                           map_table, timer);
            #endif
            chunk0_.copy_from(chunks_in.chunk1_);
            chunk1_.move_from(chunk_tmp_);
          } else {
            valid = false;
          }
        } else if (chunks_in.chunk1_ <= chunk0_) {
          if (chunks_in.chunk0_ <= chunks_in.chunk1_) {
            // [new0][new1][old0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk0_, chunks_in.chunk1_, chunk0_, connect_periodic, 2,
                           chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk0_, chunks_in.chunk1_, chunk0_, connect_periodic, 2,
                           chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
          } else if (chunks_in.chunk0_ < chunks_in.chunk1_) {
            // [new0] ... [new1][old0]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk1_, chunk0_, connect_periodic, 1, chunk_tmp_, flip_table,
                           map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk1_, chunk0_, connect_periodic, 1, chunk_tmp_, flip_table,
                           map_table, timer);
            #endif
            chunk0_.copy_from(chunks_in.chunk0_);
            chunk1_.move_from(chunk_tmp_);
          } else {
            valid = false;
          }
        } else {
          valid = false;
        }
        timer.stop(43);
      }
    } else {
      if (chunks_in.size() == 1) {
        // 2 + 1
        timer.start(44);
        const bool in0 = (orig_ <= chunk0_.range().second);
        if (chunk0_ <= chunks_in.chunk0_) {
          if (chunks_in.chunk0_ <= chunk1_) {
            // [old0][new0][old1]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunk0_, chunks_in.chunk0_, chunk1_, connect_periodic,
                           (in0  ? 0 : 2), chunk_tmp_, flip_table, map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunk0_, chunks_in.chunk0_, chunk1_, connect_periodic,
                           (in0  ? 0 : 2), chunk_tmp_, flip_table, map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
            chunk1_.clear();
          } else if (chunks_in.chunk0_ < chunk1_) {
            // [old0][new0] ... [old1]
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunk0_, chunks_in.chunk0_, 0, (in0 ? 0 : -1), chunk_tmp_, flip_table,
                           map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunk0_, chunks_in.chunk0_, 0, (in0 ? 0 : -1), chunk_tmp_, flip_table,
                           map_table, timer);
            #endif
            chunk0_.move_from(chunk_tmp_);
            if (!in0) chunk1_.fill_flip(flip_table, map_table);
          } else {
            valid = false;
          }
        } else if (chunk0_ < chunks_in.chunk0_) {
          if (chunks_in.chunk0_ <= chunk1_) {
            // [old0] ... [new0][old1]
            if (in0) chunk0_.fill_flip(flip_table, map_table);
            #ifdef LOOPER_OPENMP
            chunk_t::unify(chunks_in.chunk0_, chunk1_, 0, (in0 ? -1 : 1), chunk_tmp_, flip_table,
                           map_table, timer, estimates_tg);
            #else
            chunk_t::unify(chunks_in.chunk0_, chunk1_, 0, (in0 ? -1 : 1), chunk_tmp_, flip_table,
                           map_table, timer);
            #endif
            chunk1_.move_from(chunk_tmp_);
          } else {
            valid = false;
          }
        } else {
          valid = false;
        }
        timer.stop(44);
      } else {
        // 2 + 2
        timer.start(45);
        const bool in0 = (orig_ <= chunk0_.range().second);
        if (chunk0_ <= chunks_in.chunk0_) {
          if (chunks_in.chunk1_ <= chunk1_) {
            if (chunks_in.chunk0_ <= chunks_in.chunk1_) {
              // [old0][new0][new1][old1]
              #ifdef LOOPER_OPENMP
              chunk_t::unify(chunk0_, chunks_in.chunk0_, chunks_in.chunk1_, chunk1_,
                             connect_periodic, (in0  ? 0 : 3), chunk_tmp_, flip_table, map_table,
                             timer, estimates_tg);
              #else
              chunk_t::unify(chunk0_, chunks_in.chunk0_, chunks_in.chunk1_, chunk1_,
                             connect_periodic, (in0  ? 0 : 3), chunk_tmp_, flip_table, map_table,
                             timer);
              #endif
              chunk0_.move_from(chunk_tmp_);
              chunk1_.clear();
            } else if (chunks_in.chunk0_ < chunks_in.chunk1_) {
              // [old0][new0] ... [new1][old1]
              #ifdef LOOPER_OPENMP
              chunk_t::unify(chunk0_, chunks_in.chunk0_, connect_periodic, (in0  ? 0 : -1),
                             chunk_tmp_, flip_table, map_table, timer, estimates_tg);
              #else
              chunk_t::unify(chunk0_, chunks_in.chunk0_, connect_periodic, (in0  ? 0 : -1),
                             chunk_tmp_, flip_table, map_table, timer);
              #endif
              chunk0_.move_from(chunk_tmp_);
              #ifdef LOOPER_OPENMP
              chunk_t::unify(chunks_in.chunk1_, chunk1_, connect_periodic, (in0  ? -1 : 1),
                             chunk_tmp_, flip_table, map_table, timer, estimates_tg);
              #else
              chunk_t::unify(chunks_in.chunk1_, chunk1_, connect_periodic, (in0  ? -1 : 1),
                             chunk_tmp_, flip_table, map_table, timer);
              #endif
              chunk1_.move_from(chunk_tmp_);
            } else {
              valid = false;
            }
          } else {
            valid = false;
          }
        } else if (chunk0_ <= chunks_in.chunk1_) {
          if (chunks_in.chunk0_ <= chunk1_) {
            if (chunks_in.chunk1_ <= chunks_in.chunk0_) {
              // [old0][new1][new0][old1]
              #ifdef LOOPER_OPENMP
              chunk_t::unify(chunk0_, chunks_in.chunk1_, chunks_in.chunk0_, chunk1_,
                             connect_periodic, (in0  ? 0 : 3), chunk_tmp_, flip_table, map_table,
                             timer, estimates_tg);
              #else
              chunk_t::unify(chunk0_, chunks_in.chunk1_, chunks_in.chunk0_, chunk1_,
                             connect_periodic, (in0  ? 0 : 3), chunk_tmp_, flip_table, map_table,
                             timer);
              #endif
              chunk0_.move_from(chunk_tmp_);
              chunk1_.clear();
            } else if (chunks_in.chunk1_ < chunks_in.chunk0_) {
              // [old0][new1] ... [new0][old1]
              #ifdef LOOPER_OPENMP
              chunk_t::unify(chunk0_, chunks_in.chunk1_,
                             connect_periodic, (in0  ? 0 : -1), chunk_tmp_, flip_table, map_table,
                             timer, estimates_tg);
              #else
              chunk_t::unify(chunk0_, chunks_in.chunk1_,
                             connect_periodic, (in0  ? 0 : -1), chunk_tmp_, flip_table, map_table,
                             timer);
              #endif
              chunk0_.move_from(chunk_tmp_);
              #ifdef LOOPER_OPENMP
              chunk_t::unify(chunks_in.chunk1_, chunk0_,
                             connect_periodic, (in0  ? -1 : 1), chunk_tmp_, flip_table, map_table,
                             timer, estimates_tg);
              #else
              chunk_t::unify(chunks_in.chunk1_, chunk0_,
                             connect_periodic, (in0  ? -1 : 1), chunk_tmp_, flip_table, map_table,
                             timer);
              #endif
              chunk1_.move_from(chunk_tmp_);
            } else {
              valid = false;
            }
          } else {
            valid = false;
          }
        } else {
          valid = false;
        }
        timer.stop(45);
      }
    }
    if (!valid) {
      std::cerr << "invalid chunks in insert\n";
      boost::throw_exception(std::logic_error("chunks::insert"));
    }
    timer.stop(41);
  }

  void output_links(MPI_Comm comm, int np, int pid, int stage, int step, std::string title,
                    std::ostream& os = std::cerr) const {
    os.flush();
    MPI_Barrier(comm);
    for (int p = 0; p < np; ++p) {
      if (p == pid) {
        if (!chunk0_.empty())
          chunk0_.output_links(pid, stage, step, title + " (chunk0)", os);
        if (!chunk1_.empty())
          chunk1_.output_links(pid, stage, step, title + " (chunk1)", os);
        os.flush();
      }
      MPI_Barrier(comm);
    }
  }

  void output_collector(MPI_Comm comm, int np, int pid, int stage, int step, std::string title,
                        std::ostream& os = std::cerr) const {
    os.flush();
    MPI_Barrier(comm);
    for (int p = 0; p < np; ++p) {
      if (p == pid) {
        if (!chunk0_.empty())
          chunk0_.output_collector(pid, stage, step, title + " (chunk0)", os);
        if (!chunk1_.empty())
          chunk1_.output_collector(pid, stage, step, title + " (chunk1)", os);
        os.flush();
      }
      MPI_Barrier(comm);
    }
  }
#endif // COMMUNICATION_TEST

private:
#ifndef COMMUNICATION_TEST
  int orig_;
#endif // COMMUNICATION_TEST
  chunk_t chunk0_;
  chunk_t chunk1_;
  chunk_t chunk_tmp_;
};

} // end namespace parallel

template<typename ESTIMATE, typename COLLECTOR>
class parallel_cluster_unifier {
public:
  typedef ESTIMATE estimate_t;
  typedef COLLECTOR collector_t;
  typedef looper::union_find::node_noweight link_t;
  // typedef looper::union_find::node_c link_t;
  typedef parallel::chunks<link_t, collector_t, estimate_t> chunks_t;

  class flip_t {
  public:
    flip_t() : flip_(id_mask) {}
    void set_flip(int flip) { flip_ = flip; }
    operator int() const { return flip_; }
    bool is_open() const { return flip_ < 0; }
    void set_local_cid(int id) { flip_ = id ^ id_mask; }
    int local_cid() const { return flip_ ^ id_mask; }
    void output(std::ostream& os) const {
      if (is_open())
        os << local_cid();
      else
        os << (this->operator bool() ? "T" : "F");
    }
  private:
    static const int id_mask = -1;
    int flip_; // negative if cluster is still open
  };

  parallel_cluster_unifier(MPI_Comm comm) : comm_(comm) {}
  parallel_cluster_unifier(MPI_Comm comm, alps::parapack::timer_mpi& timer, int num_sites, std::string const& partition_str = "",
    bool duplex = true) : comm_(comm) {
    initialize(timer, num_sites, partition_str, duplex);
  }

  void initialize(alps::parapack::timer_mpi& timer, int num_sites, std::string const& partition_str = "", bool duplex = true) {
    using namespace boost::spirit;

    timer.start(36); // prepare chunks
    num_sites_ = num_sites;
    num_boundaries_ = 2 * num_sites;
    duplex_ = duplex;
    chunks_.init(num_sites_);
    chunks_r_.init(num_sites_);
    chunks_s_.init(num_sites_);

    int info;
    if ((info = MPI_Comm_size(comm_, &num_processes_)) != 0) {
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Comm_size")));
    }
    if ((info = MPI_Comm_rank(comm_, &process_id_)) != 0) {
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Comm_rank")));
    }
    timer.stop(36);

    // parse partition parameter
    timer.start(37);
    extents_.clear();
    if (partition_str.empty()) {
      extents_ = looper::prime_factorization(num_processes_); // automatic partition
    } else {
      if (!parse(partition_str.c_str(),
                 (
                  uint_p[push_back_a(extents_)] >>
                  *(punct_p >> uint_p[push_back_a(extents_)])
                 ),
                 space_p).full) {
        if (process_id_ == 0) std::cerr << "parse error: " << partition_str << std::endl;
        boost::throw_exception(std::invalid_argument("parallel_cluster_unifier: parse error"));
      }
    }
    int np = std::accumulate(extents_.begin(), extents_.end(), 1, std::multiplies<int>());
    if (np != num_processes_) {
      if (process_id_ == 0) std::cerr << "inconsistent partition: " << partition_str << std::endl;
      boost::throw_exception(std::runtime_error("parallel_cluster_unifier: inconsistent partition"));
    }
    timer.stop(37);

    // establish a cartesian topology communicator
    timer.start(38);
    coordinates_.resize(extents_.size());
    std::vector<int> periods(extents_.size(), 1); // periodic in all the directions
    if ((info = MPI_Cart_create(comm_, extents_.size(), &extents_[0], &periods[0], 0, &comm_cart_))
        != 0) {
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Cart_create")));
    }
    if ((info = MPI_Cart_coords(comm_cart_, process_id_, extents_.size(), &coordinates_[0])) != 0) {
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Cart_coords")));
    }

    if (process_id_ == 0) {
      std::cout << "Info: parallel: " << (duplex_ ? "duplex" : "simplex")
                << " mode, number of processes = " << num_processes_ << " [";
      for (int d = 0; d < extents_.size(); ++d)
        std::cout << (d != 0 ? ":" : "") << extents_[d];
      std::cout << "]\n";
    }
    timer.stop(38);

#ifndef COMMUNICATION_TEST
    // working arrays
    timer.start(39);
    flip_stage_.resize(num_boundaries_);
    timer.stop(39);
#endif

#ifdef COMMUNICATION_TEST
    int timer_id = 40;
    for (int stage = 0; stage < extents_.size(); ++stage) {  // stage loop
      const int dim = extents_.size() - stage - 1; // for transpose
      const int num_steps = (duplex_ ? extents_[dim] / 2 : extents_[dim] - 1);
      for (int step = 0; step < num_steps; ++step) {  // step loop
        timer.registrate(timer_id, "   unify:stage_" + boost::lexical_cast<std::string>(stage)
                         + ":step_" + boost::lexical_cast<std::string>(step));
        ++timer_id;
      }
    }
#endif
  }

  static void init_timer(alps::parapack::timer_mpi& timer) {
    timer.registrate(21, "  unify:all");
    timer.registrate(22, "   unify:coll.num_open_clusters");
    timer.registrate(23, "   unify:initialize_flip_table");
    timer.registrate(24, "   unify:generate_a_chunk");
    timer.registrate(25, "   unify:for_stage");
    timer.registrate(26, "    unify:chunks_r_.copy_from");
    timer.registrate(27, "    unify:for_step");
    timer.registrate(28, "     unify:chunks_s_.move_from", timer.detailed);
    timer.registrate(29, "     unify:chunks_t::sendrecv2_2way", timer.detailed | timer.start_barrier);
    timer.registrate(30, "     unify:chunks_t::sendrecv2_1way", timer.detailed | timer.start_barrier);
    timer.registrate(31, "     unify:chunks_t::sendrecv1_1way", timer.detailed | timer.start_barrier);
    timer.registrate(32, "     unify:chunks_insert", timer.detailed);
    timer.registrate(33, "     unify:update_flip_table");
    timer.registrate(34, "   unify:after_final_stage");
    timer.registrate(35, "   unify:MPI_Barrier", timer.detailed);
    timer.registrate(36, "  unifier:initialize:prepare_chunks");
    timer.registrate(37, "  unifier:initialize:partition");
    timer.registrate(38, "  unifier:initialize:cartesian");
    timer.registrate(39, "  unifier:initialize:working_arrays");
    chunks_t::init_timer(timer);
  }

  #ifdef LOOPER_OPENMP
  template<typename FRAGMENT>
  void unify(collector_t& coll, std::vector<FRAGMENT> const& fragments,
    int lower_offset, int upper_offset, std::vector<std::vector<estimate_t> >& estimates_g,
    std::vector<int>& global_id, alps::parapack::timer_mpi& timer) {
    this->unify(coll, fragments, lower_offset, upper_offset, estimates_g, timer, &global_id);
  }
  #else
  template<typename FRAGMENT>
  void unify(collector_t& coll, std::vector<FRAGMENT> const& fragments,
    int lower_offset, int upper_offset, std::vector<estimate_t>& estimatesC,
    std::vector<int>& global_id, alps::parapack::timer_mpi& timer) {
    this->unify(coll, fragments, lower_offset, upper_offset, estimatesC, timer, &global_id);
  }
  #endif

  #ifdef LOOPER_OPENMP
  template<typename FRAGMENT>
  void unify(collector_t& coll, std::vector<FRAGMENT> const& fragments,
    int lower_offset, int upper_offset, std::vector<std::vector<estimate_t> >& estimates_g,
    alps::parapack::timer_mpi& timer, std::vector<int>* global_id_ptr = 0) {
  #else
  template<typename FRAGMENT>
  void unify(collector_t& coll, std::vector<FRAGMENT> const& fragments,
    int lower_offset, int upper_offset, std::vector<estimate_t>& estimatesC,
    alps::parapack::timer_mpi& timer, std::vector<int>* global_id_ptr = 0) {
  #endif

    timer.start(21);
#ifndef COMMUNICATION_TEST
    #ifdef LOOPER_OPENMP
    std::vector<estimate_t>& estimatesC = estimates_g[0];
    #endif
#endif

    // initialize flip table
#ifndef COMMUNICATION_TEST
    timer.start(22);
    timer.stop(22);
    timer.start(23);
    const int nc_init = static_cast<int>(coll.num_clusters() + coll.num_open_clusters());
    const int noc_init = coll.num_open_clusters();
#ifdef LOOPER_OPENMP
    looper::expand(flip_, nc_init);
#else
    looper::expand(flip_, noc_init);
#endif
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static) nowait
      #endif
      for (int c = 0; c < noc_init; ++c) flip_[c].set_local_cid(c);
#ifdef LOOPER_OPENMP
      // save to_flip for [noc_init...nc_init) to flip_
      #pragma omp for schedule(static) nowait
      for (int c = noc_init; c < nc_init; ++c) flip_[c].set_flip(estimatesC[c].to_flip);
#endif
    }
    int global_id_offset = 0;
    if (global_id_ptr) {
      global_id_ptr->resize(noc_init);
      looper::expand(map_stage_, num_boundaries_);
    }
    timer.stop(23);

    // generate a chunk
    timer.start(24);
    chunks_.set(process_id_, fragments, lower_offset, upper_offset, coll, estimatesC);
    timer.stop(24);

#ifdef DEBUG_OUTPUT
    chunks_.output_links(comm_, num_processes_, process_id_, -1, -1, "initial links");
    chunks_.output_collector(comm_, num_processes_, process_id_, -1, -1, "initial collector");
#endif
#else
    int timer_id = 40;
#endif

    timer.start(25);
    for (int stage = 0; stage < extents_.size(); ++stage) {  // stage loop
#ifdef COMMUNICATION_DEBUG_OUTPUT
      if (process_id_ == 0) std::clog << "stage: " << stage << std::endl;
#endif
      // find nearest neighbor process on Cartesian coordinates
      // attention: extents_ and coordinates_ should be transposed
      const int dim = extents_.size() - stage - 1; // for transpose
      int source, dest;
      MPI_Cart_shift(comm_cart_, dim, -1, &source, &dest);
#ifdef DEBUG_OUTPUT
      output_nearest(source, dest);
#endif

      // initialize step loop
#ifndef COMMUNICATION_TEST
      timer.start(26);
      chunks_r_.copy_from(chunks_);
      timer.stop(26);
#endif

      timer.start(27);
      const int num_steps = (duplex_ ? extents_[dim] / 2 : extents_[dim] - 1);
      for (int step = 0; step < num_steps; ++step) {  // step loop
#ifdef COMMUNICATION_DEBUG_OUTPUT
        if (process_id_ == 0) std::clog << "step: " << step << std::endl;
#endif
#ifndef COMMUNICATION_TEST
        timer.start(28);
        chunks_s_.move_from(chunks_r_);
        timer.stop(28);
#endif // COMMUNICATION_TEST
#ifdef COMMUNICATION_TEST
        MPI_Barrier(comm_);
        timer.start(timer_id);
#endif
        if (duplex_) {
          if ((step + 1 < extents_[dim] / 2) || (extents_[dim] % 2 == 1)) {
            timer.start(29);
            chunks_t::sendrecv2(chunks_s_, dest, chunks_r_, source, 1, comm_cart_);
            timer.stop(29);
          } else {
            timer.start(30);
            chunks_t::sendrecv(chunks_s_, dest, chunks_r_, source, 1, comm_cart_);
            timer.stop(30);
          }
        } else {
          timer.start(31);
          chunks_t::sendrecv(chunks_s_, dest, chunks_r_, source, 1, comm_cart_);
          timer.stop(31);
        }
#ifdef COMMUNICATION_TEST
        timer.stop(timer_id);
        ++timer_id;
#endif

#ifdef DEBUG_OUTPUT
        chunks_r_.output_links(comm_, num_processes_, process_id_, stage, step, "received links");
        chunks_r_.output_collector(comm_, num_processes_, process_id_, stage, step, "received collector");
#endif

#ifndef COMMUNICATION_TEST
        bool periodic = ((stage + 1 == extents_.size()) && (step + 1 == num_steps));
        timer.start(32);
        #ifdef LOOPER_OPENMP
        chunks_.insert(chunks_r_, periodic, flip_stage_, map_stage_, timer, estimates_g);
        #else
        chunks_.insert(chunks_r_, periodic, flip_stage_, map_stage_, timer);
        #endif
        timer.stop(32);

#ifdef DEBUG_OUTPUT
        chunks_.output_links(comm_, num_processes_, process_id_, stage, step, "after unify links");
        chunks_.output_collector(comm_, num_processes_, process_id_, stage, step, "after unify collector");
        output_flips(stage, step, "flip_stage_", flip_stage_);
        output_flips(stage, step, "before update flip", flip_);
#endif

        // update flip table
        timer.start(33);
        if (global_id_ptr) {
          int nc = static_cast<int>(chunks_.get_collector().num_clusters() + chunks_.get_collector().num_open_clusters());
          int noc = chunks_.get_collector().num_open_clusters();
          int offset = global_id_offset - noc;
          global_id_offset += (nc - noc);
          #ifdef LOOPER_OPENMP
          #pragma omp parallel for schedule(static)
          #endif
          for (int c = 0; c < noc_init; ++c) {
            if (flip_[c].is_open()) {
              int cid = flip_[c].local_cid();
              flip_[c] = flip_stage_[cid];
              if (!flip_[c].is_open()) (*global_id_ptr)[c] = offset + map_stage_[cid];
            }
          }
        } else {
          #ifdef LOOPER_OPENMP
          #pragma omp parallel for schedule(static)
          #endif
          for (int c = 0; c < noc_init; ++c) {
            if (flip_[c].is_open()) {
              flip_[c] = flip_stage_[flip_[c].local_cid()];
            }
          }
        }
        timer.stop(33);

#ifdef DEBUG_OUTPUT
        output_flips(stage, step, "after update flip", flip_);
#endif
#endif // COMMUNICATION_TEST

      } // end of step loop
      timer.stop(27);

    } // end of stage loop
    timer.stop(25);

    // after final stage
#ifndef COMMUNICATION_TEST
    timer.start(34);
    // write back to to_flip
#ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    for (int c = 0; c < nc_init; ++c) estimatesC[c].to_flip = flip_[c];
#else
    for (int c = 0; c < noc_init; ++c) estimatesC[c].to_flip = flip_[c];
#endif
    coll = chunks_.get_collector();
    timer.stop(34);
    timer.start(35);
    MPI_Barrier(comm_);
    timer.stop(35);
#endif
    timer.stop(21);
  }

  template<typename T>
  void output_flips(int stage, int step, std::string title, std::vector<T> const& flips,
                    std::ostream& os = std::cerr) const;
  void output_nearest(int source, int dest, std::ostream& os = std::cerr) const;

private:
  BOOST_STATIC_ASSERT(sizeof(flip_t) == sizeof(int));
  BOOST_STATIC_ASSERT(sizeof(link_t) == sizeof(int));

  MPI_Comm comm_;
  int num_processes_;
  int process_id_;
  MPI_Comm comm_cart_; // cartesian topology communicator
  std::vector<int> extents_; // processes partition
  std::vector<int> coordinates_; // coordinates of process in Cartesian system

  int num_sites_; // number of (virtual) sites
  int num_boundaries_; // number of sites at imaginary-boundaries (2 * num_sites_)

  bool duplex_; // use duplex communication or simplex one

  // working areas
  chunks_t chunks_;
  chunks_t chunks_r_;
  chunks_t chunks_s_;

  std::vector<flip_t> flip_;
  std::vector<flip_t> flip_stage_;
  std::vector<int> map_stage_;
};


//
// implementation: debug output
//

template<typename E, typename C>
template<typename T>
void parallel_cluster_unifier<E, C>::output_flips(int stage, int step, std::string title,
  std::vector<T> const& flips, std::ostream& os) const {
  os.flush();
  MPI_Barrier(comm_);
  if (stage >= 0) {
    if (step >= 0) {
      title = std::string("stage ") + boost::lexical_cast<std::string>(stage) + ", step " +
        boost::lexical_cast<std::string>(step) + ": " + title;
    } else {
      title = std::string("stage ") + boost::lexical_cast<std::string>(stage) + ": " + title;
    }
  }
  for (int p = 0; p < num_processes_; ++p) {
    if (p == process_id_) {
      os << "process " << process_id_ << ": " << title << " [";
      for (int i = 0; i < flips.size(); ++i) {
        if (i != 0) os << ',';
        flips[i].output(os);
      }
      os << "]" << std::endl;
    }
    os.flush();
    MPI_Barrier(comm_);
  }
}

template<typename E, typename C>
void parallel_cluster_unifier<E, C>::output_nearest(int source, int dest, std::ostream& os) const {
  os.flush();
  MPI_Barrier(comm_);
  for (int p = 0; p < num_processes_; ++p) {
    if (p == process_id_) {
      os << "process " << process_id_ << ": source = " << source << ", dest = " << dest
         << std::endl;
    }
    os.flush();
    MPI_Barrier(comm_);
  }
}

} // end namespace looper

#endif // LOOPER_PARALLEL_H
