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
#include <parapack/process_mpi.h>
#include <boost/static_assert.hpp>

namespace looper {

class parallel_cluster_unifier {
public:
  class flip_t {
  public:
    flip_t() : flip_(id_mask) {}
    void set_flip(bool flip) { flip_ = (flip ? 1 : 0); }
    bool flip() const { return flip_ == 1; }
    bool open() const { return flip_ < 0; }
    void set_id(int id) { flip_ = id ^ id_mask; }
    int id() const { return flip_ ^ id_mask; }
  private:
    static const int id_mask = -1;
    int flip_; // negative if cluster is still open
  };

  parallel_cluster_unifier(alps::communicator_helper const& comm, int num_sites)
    : comm_(comm), num_sites_(num_sites), num_boundaries_(2 * num_sites),
      flip_(num_boundaries_), flip_stage_(num_boundaries_),
      links_(2 * num_boundaries_), linksD_(num_boundaries_), linksU_(num_boundaries_) {
    num_stages_ = 0;
    for (int t = (num_processes(comm_) - 1); t > 0; t = (t >> 1)) ++num_stages_;
  }

  template<typename FRAGMENT, typename RNG>
  std::vector<int> const& unify(int num_clusters, std::vector<FRAGMENT> const& fragments,
    RNG& r_uniform) const {

    // initialize local tables
    for (int c = 0; c < num_clusters; ++c) flip_[c].set_id(c);
    for (int v = 0; v < num_boundaries_; ++v) {
      if (fragments[v].is_root())
        linksD_[v].set_id(fragments[v].id()); // id = [0 ... num_clusters-1]
      else
        linksD_[v].set_parent(root_index(fragments, v));
    }

    MPI_Status status;
    for (int stage = 0; stage < num_stages_; ++stage) {
      const int stage_mask = (1 << (stage + 1)) - 1; // 000011111 for stage = 4
      const int target_mask = (1 << stage);          // 000010000

      if ((process_id(comm_) & stage_mask) == 0 &&
          process_id(comm_) + target_mask < num_processes(comm_)) {
        // master process
        const int slave = process_id(comm_) + target_mask;

        // construct link table
        for (int v = 0; v < num_boundaries_; ++v) {
          if (linksD_[v].is_root()) {
            links_[v].set_as_root();
          } else {
            links_[v].set_parent(linksD_[v].parent());
          }
        }
        MPI_Recv(&linksU_[0], num_boundaries_, MPI_INT, slave, stage, comm_.comm, &status);
        for (int v = 0; v < num_boundaries_; ++v) {
          if (linksU_[v].is_root()) {
            links_[num_boundaries_ + v].set_as_root();
          } else {
            links_[num_boundaries_ + v].set_parent(num_boundaries_ + linksU_[v].parent());
          }
        }

        // connect between lower and upper
        for (int v = 0; v < num_sites_; ++v) unify(links_, num_sites_ + v, num_boundaries_ + v);

        // at final stage bottom and top must be connected too
        if (stage + 1 == num_stages_)
          for (int v = 0; v < num_sites_; ++v) unify(links_, v, num_boundaries_ + num_sites_ + v);

        // assign cluster ID
        int nc = 0; // total number of clusters
        for (int v = 0; v < num_sites_; ++v) {
          if (links_[v].is_root) {
            links_[v].set_id(nc++);
          } else {
            int r = root_index(links_, v);
            if (r > v) {
              links_[v].set_as_root();
              links_[v].set_id(nc++);
              links_[r].set_parent(v);
            }
          }
        }
        for (int v = num_boundaries_ + num_sites_; v < 2 * num_boundaries_; ++v) {
          if (links_[v].is_root) {
            links_[v].set_id(nc++);
          } else {
            int r = root_index(links_, v);
            if ((r >= num_sites_ && r < num_boundaries_ + num_sites_) || r > v) {
              links_[v].set_as_root();
              links_[v].set_id(nc++);
              links_[r].set_parent(v);
            }
          }
        }
        int noc = nc; // number of open clusters
        for (int v = num_sites_; v < num_boundaries_ + num_sites_; ++v)
          if (links_[v].is_root) links_[v].set_id(nc++);

        // determine whether close clusters are flipped or not
        for (int c = noc; c < nc; ++c)
          flip_close_[c] = (r_uniform() < 0.5);

        // construct and send flip table for upper part
        for (int v = 0; v < num_boundaries_; ++v) {
          if (linksU_[v].is_root()) {
            const int new_id = cluster_id(links_, num_boundaries_ + v);
            if (new_id < noc)
              flip_stage_[linksU_[v].id()].set_id(new_id);
            else
              flip_stage_[linksU_[v].id()].set_flip(flip_close_[new_id]);
          }
        }
        MPI_Send(&flip_stage_[0], num_boundaries_, MPI_INT, slave, stage, comm_.comm);

        // construct flip table for lower part
        for (int v = 0; v < num_boundaries_; ++v) {
          if (linksD_[v].is_root()) {
            const int new_id = cluster_id(links_, v);
            if (new_id < noc)
              flip_stage_[linksD_[v].id()].set_id(new_id);
            else
              flip_stage_[linksD_[v].id()].set_flip(flip_close_[new_id]);
          }
        }
      } else if ((process_id(comm_) & stage_mask) == target_mask) {
        // slave process
        const int master = process_id(comm_) - target_mask;
        MPI_Send(&linksD_[0], num_boundaries_, MPI_INT, master, stage, comm_.comm);
        MPI_Recv(&flip_stage_[0], num_boundaries_, MPI_INT, master, stage, comm_.comm, &status);
      } else {
        // nothing to do
      }

      // distribute flip table to decendants
      for (int s = stage - 1; s >= 0; --s) {
        const int sm = (1 << (s + 1)) - 1;
        const int tm = (1 << s);
        if ((process_id(comm_) & sm) == 0 &&
            process_id(comm_) + tm < num_processes(comm_)) {
          const int slave = process_id(comm_) + tm;
          MPI_Send(&flip_stage_[0], num_boundaries_, MPI_INT, slave, s , comm_.comm);
        } else if ((process_id(comm_) & sm) == tm) {
          const int master = process_id(comm_) - tm;
          MPI_Recv(&flip_stage_[0], num_boundaries_, MPI_INT, master, s, comm_.comm, &status);
        } else {
          // nothing to do
        }
      }
      // update flip table
      for (int c = 0; c < num_clusters; ++c)
        if (flip_[c].open()) flip_[c] = flip_stage_[flip_[c].id()];
    }

    return flip_;
  }

private:
  BOOST_STATIC_ASSERT(sizeof(flip_t) == sizeof(int));
  BOOST_STATIC_ASSERT(sizeof(looper::union_find::node_noweight) == sizeof(int));

  const alps::communicator_helper comm_;
  const int num_sites_;                  // number of (virtual) sites
  const int num_boundaries_;             // 2 * num_sites_
  int num_stages_;                       // total number of stages

  mutable std::vector<flip_t> flip_;

  // working vectors
  mutable std::vector<bool> flip_close_;
  mutable std::vector<flip_t> flip_stage_;
  mutable std::vector<looper::union_find::node_noweight> links_;
  mutable std::vector<looper::union_find::node_noweight> linksD_;
  mutable std::vector<looper::union_find::node_noweight> linksU_;
};

} // end namespace looper

#endif // LOOPER_PARALLEL_H
