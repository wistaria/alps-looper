/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>,
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

#ifndef LOOPER_PARALLEL_H
#define LOOPER_PARALLEL_H

#include "union_find.h"
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <mpi.h>

namespace looper {

template<typename ESTIMATE, typename COLLECTOR>
class parallel_cluster_unifier {
public:
  typedef ESTIMATE estimate_t;
  typedef COLLECTOR collector_t;
  typedef looper::union_find::node_noweight link_t;

  class flip_t {
  public:
    flip_t() : flip_(id_mask) {}
    void set_flip(int flip) { flip_ = flip; }
    operator bool() const { return (flip_ & 1) == 1; }
    bool is_open() const { return flip_ < 0; }
    void set_local_cid(int id) { flip_ = id ^ id_mask; }
    int local_cid() const { return flip_ ^ id_mask; }
  private:
    static const int id_mask = -1;
    int flip_; // negative if cluster is still open
  };

  parallel_cluster_unifier(MPI_Comm comm, int num_sites) :
    comm_(comm), num_sites_(num_sites), num_boundaries_(2 * num_sites) {

    int info;
    if ((info = MPI_Comm_size(comm_, &num_processes_)) != 0) {
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Comm_size")));
    }
    if ((info = MPI_Comm_rank(comm, &process_id_)) != 0) {
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Comm_rank")));
    }
    num_stages_ = 0;
    for (int t = (num_processes_ - 1); t > 0; t = (t >> 1)) ++num_stages_;

    if (num_processes_ == 1) {
      links_.resize(num_boundaries_);
      linksD_.resize(num_boundaries_);
    } else {
      flip_stage_.resize(num_boundaries_);
      links_.resize(2 * num_boundaries_);
      linksD_.resize(num_boundaries_);
      linksU_.resize(num_boundaries_);
    }
  }

  template<typename FRAGMENT, typename RNG>
  void unify(collector_t& coll, std::vector<flip_t>& flip, std::vector<FRAGMENT> const& fragments,
    std::vector<estimate_t> const& estimates, RNG& r_uniform) {

    // initialize local tables
    const int noc_init = coll.num_open_clusters();
    for (int c = 0; c < noc_init; ++c) flip[c].set_local_cid(c);
    if (num_processes_ > 1)
      for (int c = 0; c < num_boundaries_; ++c) flip_stage_[c].set_local_cid(c);
    for (int v = 0; v < num_boundaries_; ++v) {
      if (fragments[v].is_root())
        linksD_[v].set_id(fragments[v].id()); // id = [0 ... num_open_clusters)
      else
        linksD_[v].set_parent(root_index(fragments, v));
    }
    estimatesD_.resize(coll.num_open_clusters());
    std::copy(estimates.begin(), estimates.begin() + noc_init, estimatesD_.begin());

    if (num_processes_ == 1) {

      // connect bottom and top boundaries
      std::copy(linksD_.begin(), linksD_.end(), links_.begin());
      for (int v = 0; v < num_sites_; ++v) looper::union_find::unify(links_, v, num_sites_ + v);

      // assign cluster ID
      int nc = 0;
      for (int v = 0; v < num_boundaries_; ++v) if (links_[v].is_root()) links_[v].set_id(nc++);

      estimates_.resize(nc);
      std::fill(estimates_.begin(), estimates_.end(), estimate_t());
      for (int v = 0; v < num_boundaries_; ++v)
        if (linksD_[v].is_root())
          estimates_[cluster_id(links_, v)] += estimatesD_[linksD_[v].id()];

      // accumulate cluster properties
      BOOST_FOREACH(estimate_t const& est, estimates_) coll += est;
      coll.inc_num_clusters(nc);

      // determine whether clusters are flipped or not
      flip_close_.resize(nc);
      for (int c = 0; c < nc; ++c) flip_close_[c].set_flip(r_uniform() < 0.5);
      for (int v = 0; v < num_boundaries_; ++v)
        if (linksD_[v].is_root())
          flip[linksD_[v].id()].set_flip(flip_close_[cluster_id(links_, v)]);

    } else {

      MPI_Status status;
      for (int stage = 0; stage < num_stages_; ++stage) {
        const int stage_mask = (1 << (stage + 1)) - 1; // 000011111 for stage = 4
        const int target_mask = (1 << stage);          // 000010000

        if ((process_id_ & stage_mask) == 0 && process_id_ + target_mask < num_processes_) {
          // master process
          const int slave = process_id_ + target_mask;

          std::copy(linksD_.begin(), linksD_.end(), links_.begin());
          MPI_Recv(&coll_buf_, sizeof(collector_t), MPI_BYTE, slave, stage, comm_, &status);
          coll += coll_buf_;
          estimatesU_.resize(coll_buf_.num_open_clusters());
          MPI_Recv(&linksU_[0], num_boundaries_, MPI_INT, slave, stage, comm_, &status);
          MPI_Recv(&estimatesU_[0], sizeof(estimate_t) * coll_buf_.num_open_clusters(), MPI_BYTE,
                   slave, stage, comm_, &status);
          for (int v = 0; v < num_boundaries_; ++v) {
            if (linksU_[v].is_root()) {
              links_[num_boundaries_ + v] = link_t();
            } else {
              links_[num_boundaries_ + v].set_parent(num_boundaries_ + linksU_[v].parent());
            }
          }

          // connect between lower and upper
          for (int v = 0; v < num_sites_; ++v)
            looper::union_find::unify(links_, num_sites_ + v, num_boundaries_ + v);

          // at final stage bottom and top boundaries must be connected too
          if (stage + 1 == num_stages_) {
            for (int v = 0; v < num_sites_; ++v)
              looper::union_find::unify(links_, v, num_boundaries_ + num_sites_ + v);
          } else {
            for (int v = 0; v < num_sites_; ++v) set_root(links_, v);
            for (int v = num_boundaries_ + num_sites_; v < 2 * num_boundaries_; ++v)
              set_root(links_, v);
          }

          // assign cluster ID
          int nc = 0; // total number of clusters
          for (int v = 0; v < num_sites_; ++v)
            if (links_[v].is_root()) links_[v].set_id(nc++);
          for (int v = num_boundaries_ + num_sites_; v < 2 * num_boundaries_; ++v)
            if (links_[v].is_root()) links_[v].set_id(nc++);
          int noc = nc; // number of open clusters
          for (int v = num_sites_; v < num_boundaries_ + num_sites_; ++v)
            if (links_[v].is_root()) links_[v].set_id(nc++);
          if (stage + 1 == num_stages_) noc = 0; // at finil stage there's no open clusters

          // determine whether close clusters are flipped or not
          flip_close_.resize(nc);
          for (int c = noc; c < nc; ++c) flip_close_[c].set_flip(r_uniform() < 0.5);
          estimates_.resize(nc);
          std::fill(estimates_.begin(), estimates_.end(), estimate_t());
          for (int v = 0; v < num_boundaries_; ++v) {
            if (linksD_[v].is_root())
              estimates_[cluster_id(links_, v)] += estimatesD_[linksD_[v].id()];
            if (linksU_[v].is_root())
              estimates_[cluster_id(links_, num_boundaries_ + v)] += estimatesU_[linksU_[v].id()];
          }

          // construct and send flip table for upper part
          for (int v = 0; v < num_boundaries_; ++v) {
            if (linksU_[v].is_root()) {
              const int new_id = cluster_id(links_, num_boundaries_ + v);
              const int old_id = linksU_[v].id();
              if (new_id < noc)
                flip_stage_[old_id].set_local_cid(new_id);
              else
                flip_stage_[old_id].set_flip(flip_close_[new_id]);
            }
          }
          MPI_Send(&flip_stage_[0], num_boundaries_, MPI_INT, slave, stage, comm_);

          // construct flip table for lower part
          for (int v = 0; v < num_boundaries_; ++v) {
            if (linksD_[v].is_root()) {
              const int new_id = cluster_id(links_, v);
              const int old_id = linksD_[v].id();
              if (new_id < noc)
                flip_stage_[old_id].set_local_cid(new_id);
              else
                flip_stage_[old_id].set_flip(flip_close_[new_id]);
            }
          }

          // distribute flip table to decendants
          distribute(stage);

          // accumulate measurements of closed clusters
          for (int c = noc; c < nc; ++c) coll += estimates_[c];
          coll.set_num_open_clusters(noc);
          coll.inc_num_clusters(nc - noc);

          // update links and estimates (linksD_ and estimatesD_)
          if (stage + 1 != num_stages_) {
            for (int vd = 0; vd < num_boundaries_; ++vd) {
              const int v = (vd < num_sites_ ? vd : num_boundaries_ + vd);
              if (links_[v].is_root()) {
                linksD_[vd].set_id(links_[v].id());
              } else {
                int r = root_index(links_, v);
                linksD_[vd].set_parent((r < num_sites_ ? r : r - num_boundaries_));
              }
            }
            estimatesD_.resize(noc);
            std::copy(estimates_.begin(), estimates_.begin() + noc, estimatesD_.begin());
          }

        } else if ((process_id_ & stage_mask) == target_mask) {
          // slave process
          const int master = process_id_ - target_mask;
          MPI_Send(&coll, sizeof(collector_t), MPI_BYTE, master, stage, comm_);
          MPI_Send(&linksD_[0], num_boundaries_, MPI_INT, master, stage, comm_);
          MPI_Send(&estimatesD_[0], coll.num_open_clusters() * sizeof(estimate_t), MPI_BYTE,
                   master, stage, comm_);
          MPI_Recv(&flip_stage_[0], num_boundaries_, MPI_INT, master, stage, comm_, &status);
          // distribute flip table to decendants
          distribute(stage);
        } else {
          // distribute flip table to decendants
          distribute(stage);
        }

        // update flip table
        for (int c = 0; c < noc_init; ++c)
          if (flip[c].is_open()) flip[c] = flip_stage_[flip[c].local_cid()];
      }
    }
  }

protected:
  // distribute flip table to decendants
  void distribute(int stage) {
    MPI_Status status;
    for (int s = stage - 1; s >= 0; --s) {
      const int sm = (1 << (s + 1)) - 1;
      const int tm = (1 << s);
      if ((process_id_ & sm) == 0 && process_id_ + tm < num_processes_) {
        const int slave = process_id_ + tm;
        MPI_Send(&flip_stage_[0], num_boundaries_, MPI_INT, slave, s , comm_);
      } else if ((process_id_ & sm) == tm) {
        const int master = process_id_ - tm;
        MPI_Recv(&flip_stage_[0], num_boundaries_, MPI_INT, master, s, comm_, &status);
      } else {
        // nothing to do
      }
    }
  }

private:
  BOOST_STATIC_ASSERT(sizeof(flip_t) == sizeof(int));
  BOOST_STATIC_ASSERT(sizeof(looper::union_find::node_noweight) == sizeof(int));

  MPI_Comm comm_;
  int num_processes_;
  int process_id_;
  int num_sites_; // number of (virtual) sites
  int num_boundaries_; // number of sites at imaginary-boundaries (2 * num_sites_)
  int num_stages_; // total number of stages

  // working areas
  collector_t coll_buf_;
  std::vector<flip_t> flip_close_;
  std::vector<flip_t> flip_stage_;
  std::vector<link_t> links_;
  std::vector<link_t> linksD_;
  std::vector<link_t> linksU_;
  std::vector<estimate_t> estimates_;
  std::vector<estimate_t> estimatesD_;
  std::vector<estimate_t> estimatesU_;
};

} // end namespace looper

#endif // LOOPER_PARALLEL_H
