/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOP_WORKER_PI_H
#define LOOP_WORKER_PI_H

#include "loop_worker.h"

struct cluster_info_pi
{
  cluster_info_pi(bool t = false)
    : to_flip(t), mag0(0), size(0), mag(0), length(0) {}
  bool to_flip;
  int mag0;
  int size;
  double mag;
  double length;
};

enum operator_type { bond_diagonal = 0, bond_offdiagonal = 1 };

struct local_operator_t {
  local_operator_t() {}
  local_operator_t(operator_type tp, int loc, double t)
    : type(tp), location(loc), loop0(), loop1(), time_(t) {}
  void flip() { type = operator_type(type ^ 1); }
  bool is_diagonal() const { return !is_offdiagonal(); }
  bool is_offdiagonal() const { return type & 1; }
  unsigned int pos() const { return location; }
  double time() const { return time_; }
  operator_type type;
  unsigned int location;
  unsigned int loop0, loop1;
  double time_;
  bool is_bond() const { return true; }
  bool is_locked() const { return true; }
};

inline alps::ODump& operator<<(alps::ODump& od, const local_operator_t& op)
{
  int t = op.type;
  od << t << op.location << op.time_;
  return od;
}

inline alps::IDump& operator>>(alps::IDump& id, local_operator_t& op)
{
  int t;
  id >> t >> op.location >> op.time_;
  op.type = operator_type(t);
  return id;
}

class qmc_worker_pi : public qmc_worker_base
{
public:
  typedef looper::path_integral                         qmc_type;
  typedef qmc_worker_base                               super_type;
  typedef loop_config::graph_type                       graph_type;
  typedef loop_config::local_graph                      local_graph;
  typedef looper::local_operator<qmc_type, local_graph> local_operator;
  // typedef local_operator_t local_operator;
  typedef looper::union_find::node                      cluster_fragment;
  typedef cluster_info_pi                               cluster_info;

  qmc_worker_pi(const alps::ProcessList& w, const alps::Parameters& p, int n);

  void dostep();

  void save(alps::ODump& dp) const {
    super_type::save(dp);
    dp << spins << operators;
  }
  void load(alps::IDump& dp) {
    super_type::load(dp);
    dp >> spins >> operators;
  }

private:
  double beta;
  looper::model_parameter mp;
  double energy_offset;

  looper::virtual_lattice<graph_type> vlat;
  bool is_bipartite;

  looper::graph_chooser<local_graph, super_type::engine_type> chooser;

  std::vector<int> spins;
  std::vector<local_operator> operators;

  // working area
  std::vector<int> spins_c;
  std::vector<local_operator> operators_p;
  std::vector<cluster_fragment> fragments;
  std::vector<int> current;
  std::vector<cluster_info> clusters;
};

#endif // LOOP_WORKER_PI_H
