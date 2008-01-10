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

#ifndef LOOPER_PARALLEL_H
#define LOOPER_PARALLEL_H

#include <looper/union_find.h>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <mpi.h>

namespace looper {

template<typename ESTIMATE, typename ACCUMULATE>
class parallel_cluster_unifier {
public:
  typedef ESTIMATE estimate_t;
  typedef ACCUMULATE accumulate_t;
  typedef looper::union_find::node_noweight link_t;

  class flip_t {
  public:
    flip_t() : flip_(id_mask) {}
    void set_flip(int flip) { flip_ = flip; }
    bool flip() const { return flip_ == 1; }
    bool open() const { return flip_ < 0; }
    void set_id(int id) { flip_ = id ^ id_mask; }
    int id() const { return flip_ ^ id_mask; }
  private:
    static const int id_mask = -1;
    int flip_; // negative if cluster is still open
  };

  struct info_t {
    accumulate_t accum;
    int nc;
  };

  parallel_cluster_unifier(MPI_Comm comm, int num_sites) :
    comm_(comm), num_sites_(num_sites), num_boundaries_(2 * num_sites) {

    int info;
    if ((info = MPI_Comm_size(comm_, &num_processes_)) != 0)
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Comm_size")));
    if ((info = MPI_Comm_rank(comm, &process_id_)) != 0)
      boost::throw_exception(std::runtime_error(("Error " +
        boost::lexical_cast<std::string>(info) + " in MPI_Comm_rank")));
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
  std::vector<flip_t> const& unify(int num_clusters, std::vector<FRAGMENT> const& fragments,
    std::vector<estimate_t> const& estimates, accumulate_t& accum, RNG& r_uniform) {

    // initialize local tables
    flip_.resize(num_clusters);
    for (int c = 0; c < num_clusters; ++c) flip_[c].set_id(c);
    for (int v = 0; v < num_boundaries_; ++v) {
      if (fragments[v].is_root())
        linksD_[v].set_id(fragments[v].id()); // id = [0 ... num_clusters-1]
      else
        linksD_[v].set_parent(root_index(fragments, v));
    }
    estimatesD_.resize(num_clusters);
    std::copy(estimates.begin(), estimates.begin() + num_clusters, estimatesD_.begin());
    info_.nc = num_clusters;
    info_.accum = accum;

    if (num_processes_ == 1) {

      std::copy(linksD_.begin(), linksD_.end(), links_.begin());
      for (int v = 0; v < num_sites_; ++v) looper::union_find::unify(links_, v, num_sites_ + v);

      // assign cluster ID
      int nc = 0;
      for (int v = 0; v < num_boundaries_; ++v) if (links_[v].is_root()) links_[v].set_id(nc++);

      estimates_.resize(nc);
      std::fill(estimates_.begin(), estimates_.end(), estimate_t());
      for (int v = 0; v < num_boundaries_; ++v)
        estimates_[cluster_id(links_, v)] += estimatesD_[cluster_id(linksD_, v)];

      // accumulate cluster properties
      BOOST_FOREACH(estimate_t const& est, estimates_) accum += est;

      // determine whether clusters are flipped or not
      BOOST_FOREACH(flip_t& f, flip_) f.set_flip(r_uniform() < 0.5);

    } else {

      MPI_Status status;
      for (int stage = 0; stage < num_stages_; ++stage) {
        const int stage_mask = (1 << (stage + 1)) - 1; // 000011111 for stage = 4
        const int target_mask = (1 << stage);          // 000010000

        if ((process_id_ & stage_mask) == 0 &&
            process_id_ + target_mask < num_processes_) {
          // master process
          const int slave = process_id_ + target_mask;

          // construct link table
          for (int v = 0; v < num_boundaries_; ++v) {
            if (linksD_[v].is_root()) {
              links_[v] = link_t();
            } else {
              links_[v].set_parent(linksD_[v].parent());
            }
          }
          if (stage != 0) {
            estimatesD_.resize(info_.nc);
            std::copy(estimates_.begin(), estimates_.begin() + info_.nc, estimatesD_.begin());
          }
            
          MPI_Recv(&info_, sizeof(info_t), MPI_BYTE, slave, stage, comm_, &status);
          accum += info_.accum;
          estimatesU_.resize(info_.nc);
          MPI_Recv(&linksU_[0], num_boundaries_, MPI_INT, slave, stage, comm_, &status);
          MPI_Recv(&estimatesU_[0], sizeof(estimate_t) * info_.nc, MPI_BYTE, slave, stage, comm_,
                   &status);
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

          // at final stage bottom and top must be connected too
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
          for (int c = noc; c < nc; ++c) flip_close_[c] = (r_uniform() < 0.5);
          estimates_.resize(nc);
          std::fill(estimates_.begin(), estimates_.end(), estimate_t());

          // construct and send flip table for upper part
          for (int v = 0; v < num_boundaries_; ++v) {
            const int new_id = cluster_id(links_, num_boundaries_ + v);
            const int old_id = cluster_id(linksU_, v);
            if (linksU_[v].is_root()) {
              if (new_id < noc)
                flip_stage_[old_id].set_id(new_id);
              else
                flip_stage_[old_id].set_flip(flip_close_[new_id]);
            }
            estimates_[new_id] += estimatesU_[old_id];
          }
          MPI_Send(&flip_stage_[0], num_boundaries_, MPI_INT, slave, stage, comm_);

          // construct flip table for lower part
          for (int v = 0; v < num_boundaries_; ++v) {
            const int new_id = cluster_id(links_, v);
            const int old_id = cluster_id(linksD_, v);
            if (linksD_[v].is_root()) {
              if (new_id < noc)
                flip_stage_[old_id].set_id(new_id);
              else
                flip_stage_[old_id].set_flip(flip_close_[new_id]);
            }
            estimates_[new_id] += estimatesD_[old_id];
          }

          // accumulate measurements of closed clusters
          for (int c = noc; c < nc; ++c) accum += estimates_[c];
          info_.nc = noc;
          info_.accum = accum;

        } else if ((process_id_ & stage_mask) == target_mask) {
          // slave process
          const int master = process_id_ - target_mask;
          MPI_Send(&info_, sizeof(info_t), MPI_BYTE, master, stage, comm_);
          MPI_Send(&linksD_[0], num_boundaries_, MPI_INT, master, stage, comm_);
          if (stage == 0) {
            MPI_Send(&estimatesD_[0], info_.nc * sizeof(estimate_t), MPI_BYTE, master, stage,
              comm_);
          } else {
            MPI_Send(&estimates_[0], info_.nc * sizeof(estimate_t), MPI_BYTE, master, stage,
              comm_);
          }
          MPI_Recv(&flip_stage_[0], num_boundaries_, MPI_INT, master, stage, comm_, &status);
        } else {
          // nothing to do
        }

        // distribute flip table to decendants
        for (int s = stage - 2; s >= 0; --s) {
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
        // update flip table
        for (int c = 0; c < num_clusters; ++c)
          if (flip_[c].open()) flip_[c] = flip_stage_[flip_[c].id()];
      }

    }

    return flip_;
  }

private:
  BOOST_STATIC_ASSERT(sizeof(flip_t) == sizeof(int));
  BOOST_STATIC_ASSERT(sizeof(looper::union_find::node_noweight) == sizeof(int));

  MPI_Comm comm_;
  int num_processes_;
  int process_id_;
  int num_sites_;      // number of (virtual) sites
  int num_boundaries_; // number of sites at imaginary-boundaries (2 * num_sites_)
  int num_stages_;     // total number of stages

  std::vector<flip_t> flip_;

  // working areas
  info_t info_;
  std::vector<int> flip_close_;
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
