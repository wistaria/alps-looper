/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_config.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/expand.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/type.h>
#include <looper/union_find.h>
#include <alps/parapack/worker.h>
#include <alps/parapack/exchange.h>
#include <alps/numeric/is_zero.hpp>

namespace {

class loop_worker : public alps::parapack::mc_worker, private loop_config {
public:
  typedef looper::path_integral mc_type;

  typedef looper::local_operator<mc_type, loop_graph_t, time_t> local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef operator_string_t::iterator operator_iterator;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef looper::estimator<measurement_set, mc_type, lattice_t, time_t>::type estimator_t;
  typedef estimator_t::improved_estimator::estimate estimate_t;
  typedef estimator_t::normal_estimator::estimate normal_estimate_t;
  typedef double weight_parameter_type;

  typedef boost::exponential_distribution<> expdist_t;

  loop_worker(alps::Parameters const& p);
  virtual ~loop_worker() {}

  void init_observables(alps::Parameters const& params, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs);

  // for exchange Monte Carlo
  void set_beta(double beta) { temperature.set_beta(beta); }
  double weight_parameter() const { return operators.size(); }
  static double log_weight(double gw, double beta) { return std::log(beta) * gw; }

  void save(alps::ODump& dp) const { dp << mcs << spins << operators; }
  void load(alps::IDump& dp) { dp >> mcs >> spins >> operators; }

protected:
  template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
  void dispatch(alps::ObservableSet& obs, COLLECTOR& coll, std::vector<ESTIMATE>& estimates);

private:
  // helpers
  lattice_t lattice;
  model_t model;

  // parameters
  looper::temperature temperature;
  double beta;
  bool enable_improved_estimator;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // observables
  estimator_t estimator;
  estimator_t::improved_estimator::collector coll_i;
  estimator_t::normal_estimator::collector coll_n;

  // working vectors
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<double> times;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<cluster_info_t> clusters;
  std::vector<estimator_t::improved_estimator::estimate> estimates_i;
  std::vector<estimator_t::normal_estimator::estimate> estimates_n;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::Parameters const& p)
  : alps::parapack::mc_worker(p), lattice(p), model(p, lattice, /* is_path_integral = */ true),
    temperature(p), mcs(p) {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(support_longitudinal_field, support_negative_sign);

  enable_improved_estimator = (!model.has_field()) && (!p.defined("DISABLE_IMPROVED_ESTIMATOR"));
  if (!enable_improved_estimator) std::cerr << "WARNING: improved estimator is disabled\n";

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);
  perm.resize(max_virtual_sites(lattice));

  // initialize estimators
  estimator.initialize(p, lattice, model.is_signed(), enable_improved_estimator);
}

void loop_worker::init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
  obs << make_observable(alps::SimpleRealObservable("Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Volume"));
  obs << make_observable(alps::SimpleRealObservable("Number of Sites"));
  obs << make_observable(alps::SimpleRealObservable("Number of Clusters"));
  if (model.is_signed()) {
    obs << alps::RealObservable("Sign");
    if (enable_improved_estimator) {
      obs << alps::RealObservable("Weight of Zero-Meron Sector");
      obs << alps::RealObservable("Sign in Zero-Meron Sector");
    }
  }
  estimator.init_observables(obs, model.is_signed(), enable_improved_estimator);
}

void loop_worker::run(alps::ObservableSet& obs) {
  // if (!mcs.can_work()) return;
  ++mcs;
  beta = 1.0 / temperature(mcs());

  //       FIELD               SIGN                IMPROVE
  dispatch<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::false_>(obs, coll_n, estimates_n);
  dispatch<boost::mpl::true_,  boost::mpl::false_, boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::true_,  boost::mpl::false_, boost::mpl::false_>(obs, coll_n, estimates_n);
  dispatch<boost::mpl::false_, boost::mpl::true_,  boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::false_, boost::mpl::true_,  boost::mpl::false_>(obs, coll_n, estimates_n);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs, coll_n, estimates_n);
}


template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
void loop_worker::dispatch(alps::ObservableSet& obs, COLLECTOR& coll,
  std::vector<ESTIMATE>& estimates) {
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      enable_improved_estimator != IMPROVE()) return;

  int nrs = num_sites(lattice.rg());
  int nvs = num_sites(lattice.vg());

  //
  // diagonal update and cluster construction
  //

  // initialize spin & operator information
  std::swap(operators, operators_p); operators.resize(0);
  // insert a diagonal operator at the end of operators_p
  operators_p.push_back(local_operator_t(0, local_operator_t::location_t(), 1));
  for (int s = 0; s < nvs; ++s) spins_c[s] = spins[s];

  // fill times
  expdist_t expdist(beta * model.graph_weight());
  times.resize(0);
  double t = 0;
  while (t < 1) {
    t += expdist(generator_01());
    times.push_back(t);
  } // a sentinel (t >= 1) will be appended

  // initialize cluster information (setup cluster fragments)
  int fragment_offset = nvs;
  looper::expand(fragments, fragment_offset + operators_p.size() + times.size());
  for (int s = 0; s < nvs; ++s) {
    fragments[s] = cluster_fragment_t();
    current[s] = s;
  }

  coll.reset(estimator);
  looper::normal_accumulator<estimator_t, IMPROVE> accum_n(coll, lattice, estimator);
  for (int s = 0; s < nvs; ++s) accum_n.start_bottom(time_t(0), s, spins_c[s]);
  int negop = 0; // number of operators with negative weights

  int fid = fragment_offset;
  std::vector<double>::iterator tmi = times.begin();
  for (operator_iterator opi = operators_p.begin(); opi != operators_p.end();) {
    // diagonal update & labeling
    if (*tmi < opi->time()) {
      loop_graph_t g = model.choose_graph(generator_01());
      if ((is_bond(g) && is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                                          spins_c[target(pos(g), lattice.vg())])) ||
          (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
        operators.push_back(local_operator_t(g, *tmi));
        ++tmi;
      } else {
        ++tmi;
        continue;
      }
    } else {
      if (opi->is_diagonal()) {
        ++opi;
        continue;
      } else {
        operators.push_back(*opi);
        ++opi;
      }
    }

    operator_iterator oi = operators.end() - 1;
    if (oi->is_bond()) {
      int b = oi->pos();
      int s0 = source(b, lattice.vg());
      int s1 = target(b, lattice.vg());
      if (oi->is_offdiagonal()) {
        oi->assign_graph(model.choose_offdiagonal(generator_01(), oi->loc(),
          spins_c[s0], spins_c[s1]));
        accum_n.end_b(oi->time(), b, s0, s1, spins_c[s0], spins_c[s1]);
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
        accum_n.begin_b(oi->time(), b, s0, s1, spins_c[s0], spins_c[s1]);
        if (SIGN()) negop += model.bond_sign(oi->pos());
      }
      boost::tie(fid, current[s0], current[s1], oi->loop0, oi->loop1) =
        reconnect(fragments, fid, oi->graph(), current[s0], current[s1]);
    } else {
      int s = oi->pos();
      if (oi->is_offdiagonal()) {
        accum_n.end_s(oi->time(), s, spins_c[s]);
        spins_c[s] ^= 1;
        accum_n.begin_s(oi->time(), s, spins_c[s]);
        if (SIGN()) negop += model.site_sign(oi->pos());
      }
      boost::tie(fid, current[s], oi->loop0, oi->loop1) =
        reconnect(fragments, fid, oi->graph(), current[s]);
    }
  }
  for (int s = 0; s < nvs; ++s) accum_n.stop_top(time_t(1), s, spins_c[s]);
  double sign = ((negop & 1) == 1) ? -1 : 1;
  int num_fragments = fid - fragment_offset;

  // symmetrize spins
  if (max_virtual_sites(lattice) == 1) {
    for (int s = 0; s < nvs; ++s) unify(fragments, s, current[s]);
  } else {
    for (int rs = 0; rs < nrs; ++rs) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      boost::tie(vsi, vsi_end) = sites(lattice, rs);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins_c.begin() + offset, generator_01());
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, current[offset+perm[i]]);
    }
  }

  //
  // cluster flip
  //

  // assign cluster id
  int nc = set_id(fragments, 0, fragment_offset + num_fragments, 0);
  copy_id(fragments, 0, fragment_offset + num_fragments);

  // accumulate physical property of clusters
  looper::expand(clusters, nc);
  looper::expand(estimates, nc);
  for (int c = 0; c < nc; ++c) {
    clusters[c] = cluster_info_t();
    estimates[c].reset(estimator);
  }
  if (IMPROVE() || FIELD()) {
    for (int s = 0; s < nvs; ++s) spins_c[s] = spins[s];
    cluster_info_t::accumulator<cluster_fragment_t, FIELD, SIGN, IMPROVE>
      weight(clusters, fragments, model.field(), model.bond_sign(), model.site_sign());
    looper::improved_accumulator<estimator_t, cluster_fragment_t, ESTIMATE, IMPROVE>
      accum_i(estimates, lattice, estimator, fragments);
    for (int s = 0; s < nvs; ++s) {
      weight.start_bottom(s, time_t(0), s, spins_c[s]);
      accum_i.start_bottom(s, time_t(0), s, spins_c[s]);
    }
    for (operator_iterator opi = operators.begin(); opi != operators.end(); ++opi) {
      time_t t = opi->time();
      if (opi->is_bond()) {
        if (!opi->is_frozen_bond_graph()) {
          int b = opi->pos();
          int s0 = source(b, lattice.vg());
          int s1 = target(b, lattice.vg());
          weight.end_b(opi->loop_l0(), opi->loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
          accum_i.end_b(opi->loop_l0(), opi->loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
          if (opi->is_offdiagonal()) {
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
          }
          weight.begin_b(opi->loop_u0(), opi->loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
          accum_i.begin_b(opi->loop_u0(), opi->loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        }
      } else {
        if (!opi->is_frozen_site_graph()) {
          int s = opi->pos();
          weight.end_s(opi->loop_l(), t, s, spins_c[s]);
          accum_i.end_s(opi->loop_l(), t, s, spins_c[s]);
          if (opi->is_offdiagonal()) spins_c[s] ^= 1;
          weight.begin_s(opi->loop_u(), t, s, spins_c[s]);
          accum_i.begin_s(opi->loop_u(), t, s, spins_c[s]);
        }
      }
    }
    for (int s = 0; s < nvs; ++s) {
      weight.stop_top(current[s], time_t(1), s, spins_c[s]);
      accum_i.stop_top(current[s], time_t(1), s, spins_c[s]);
    }
  }

  // accumulate cluster properties
  coll.set_num_operators(operators.size());
  if (IMPROVE()) for (int c = 0; c < nc; ++c) coll += estimates[c];
  coll.set_num_clusters(nc);

  // determine whether clusters are flipped or not
  double improved_sign = sign;
  for (int c = 0; c < nc; ++c) {
    estimates[c].to_flip =
      ((2*uniform_01()-1) < (FIELD() ? std::tanh(beta * clusters[c].weight) : 0));
    if (SIGN() && IMPROVE() && (clusters[c].sign & 1) == 1) improved_sign = 0;
  }

  // flip operators & spins
  for (operator_iterator opi = operators.begin(); opi != operators.end(); ++opi) {
    if (estimates[fragments[opi->loop_0()].id()].to_flip ^
        estimates[fragments[opi->loop_1()].id()].to_flip) opi->flip();
  }
  for (int s = 0; s < nvs; ++s)
    if (estimates[fragments[s].id()].to_flip) spins[s] ^= 1;

  //
  // measurement
  //

  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Volume"] << (double)lattice.volume();
  obs["Number of Sites"] << (double)nrs;
  obs["Number of Clusters"] << coll.num_clusters();
  if (SIGN()) {
    if (IMPROVE()) {
      obs["Sign"] << improved_sign;
      if (alps::numeric::is_zero(improved_sign)) {
        obs["Weight of Zero-Meron Sector"] << 0.;
      } else {
        obs["Weight of Zero-Meron Sector"] << 1.;
        obs["Sign in Zero-Meron Sector"] << improved_sign;
      }
    } else {
      obs["Sign"] << sign;
    }
  }
  double nop = coll.num_operators();
  double ene = model.energy_offset() - nop / beta;
  if (FIELD())
    for (int c = 0; c < nc; ++c)
      ene += (estimates[c].to_flip ? -clusters[c].weight : clusters[c].weight);
  coll.set_energy(ene);

  coll.commit(obs, estimator, lattice, beta, improved_sign, nop);
}

typedef looper::evaluator<loop_config::measurement_set> loop_evaluator;

//
// dynamic registration to the factories
//

PARAPACK_REGISTER_ALGORITHM(loop_worker, "loop");
PARAPACK_REGISTER_ALGORITHM(loop_worker, "loop; path integral");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<loop_worker>, "loop; exchange");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<loop_worker>, "loop; path integral; exchange");
PARAPACK_REGISTER_EVALUATOR(loop_evaluator, "loop");

} // end namespace
