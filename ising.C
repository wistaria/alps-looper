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
#include <looper/permutation.h>
#include <looper/union_find.h>
#include <looper/temperature.h>
#include <alps/parapack/worker.h>
#include <alps/parapack/exchange.h>

#ifndef LOOPER_ONLY_PATH_INTEGRAL

namespace {

struct dummy_operator {
  static bool is_offdiagonal() { return false; }
  static int pos() { return 0; }
  static bool is_site() { return false; }
  static bool is_bond() { return false; }
};

class loop_worker : public alps::parapack::mc_worker, private loop_config {
public:
  typedef looper::classical mc_type;

  typedef std::vector<dummy_operator> operator_string_t;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef looper::estimator<measurement_set, mc_type, lattice_t, time_t>::type estimator_t;
  typedef estimator_t::improved_estimator::estimate estimate_t;
  typedef estimator_t::improved_estimator::minimal_estimate minimal_estimate_t;
  typedef double weight_parameter_type;

  loop_worker(alps::Parameters const& p);
  virtual ~loop_worker() {}

  void init_observables(alps::Parameters const& params, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs);

  // for exchange Monte Carlo
  void set_beta(double beta) { temperature.set_beta(beta); }
  double weight_parameter() const { return nop; }
  static double log_weight(double gw, double beta) { return std::log(beta) * gw; }

  void save(alps::ODump& dp) const { dp << mcs << spins; }
  void load(alps::IDump& dp) { dp >> mcs >> spins; }

protected:
  template<typename FIELD, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
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

  // observables
  estimator_t estimator;
  estimator_t::improved_estimator::collector coll_i;
  estimator_t::normal_estimator::collector coll_n;

  // working vectors
  int nop;
  std::vector<cluster_fragment_t> fragments;
  std::vector<cluster_info_t> clusters;
  std::vector<estimate_t> estimates_i;
  std::vector<minimal_estimate_t> estimates_m;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::Parameters const& p)
  : alps::parapack::mc_worker(p), lattice(p), model(p, lattice), temperature(p), mcs(p) {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(support_longitudinal_field, /* support_negative_sign = */ false);
  if (model.is_quantal())
    boost::throw_exception(std::invalid_argument("not classical Ising model"));

  enable_improved_estimator = !model.has_field() && !p.defined("DISABLE_IMPROVED_ESTIMATOR");
  if (!enable_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
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
  estimator.init_observables(obs, false, enable_improved_estimator);
}

void loop_worker::run(alps::ObservableSet& obs) {
  // if (!mcs.can_work()) return;
  ++mcs;
  beta = 1.0 / temperature(mcs());

  //       FIELD               IMPROVE
  dispatch<boost::mpl::true_,  boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::true_,  boost::mpl::false_>(obs, coll_n, estimates_m);
  dispatch<boost::mpl::false_, boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::false_, boost::mpl::false_>(obs, coll_n, estimates_m);
}


template<typename FIELD, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
void loop_worker::dispatch(alps::ObservableSet& obs, COLLECTOR& coll,
  std::vector<ESTIMATE>& estimates) {
  if (model.has_field() != FIELD() ||
      enable_improved_estimator != IMPROVE()) return;

  int nrs = num_sites(lattice.rg());
  int nvs = num_sites(lattice.vg());

  //
  // cluster construction
  //

  // initialize cluster information (setup cluster fragments)
  fragments.resize(0); fragments.resize(nvs);

  // initialize measurements
  coll.reset(estimator);
  looper::normal_accumulator<estimator_t, IMPROVE> accum_n(coll, lattice, estimator);
  if (!IMPROVE()) {
    for (int s = 0; s < nvs; ++s) {
      accum_n.start_bottom(time_t(0), s, spins[s]);
      accum_n.stop_top(time_t(1), s, spins[s]);
    }
  }

  // building up clusters
  nop = 0;
  boost::variate_generator<engine_type&, boost::exponential_distribution<> >
    r_time(engine(), boost::exponential_distribution<>(beta * model.graph_weight()));
  double t = r_time();
  while (t < 1) {
    loop_graph_t g = model.choose_graph(generator_01());
    int s0 = source(pos(g), lattice.vg());
    int s1 = target(pos(g), lattice.vg());
    if (is_compatible(g, spins[s0], spins[s1])) {
      unify(fragments, s0, s1);
      ++nop;
    }
    t += r_time();
  }

  // symmetrize spins
  if (max_virtual_sites(lattice) > 1) {
    for (int rs = 0; rs < nrs; ++rs) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      boost::tie(vsi, vsi_end) = sites(lattice, rs);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins.begin() + offset, generator_01());
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, offset+perm[i]);
    }
  }

  //
  // cluster flip
  //

  // assign cluster id
  int nc = set_id(fragments, 0, fragments.size(), 0);
  copy_id(fragments, 0, fragments.size());

  // accumulate physical property of clusters
  looper::expand(clusters, nc);
  looper::expand(estimates, nc);
  for (int c = 0; c < nc; ++c) {
    clusters[c] = cluster_info_t();
    estimates[c].reset(estimator);
  }

  if (IMPROVE() || FIELD()) {
    cluster_info_t::accumulator<cluster_fragment_t, FIELD,
      boost::mpl::false_, IMPROVE> weight(clusters, fragments, model.field(), 0, 0);
    looper::improved_accumulator<estimator_t, cluster_fragment_t, ESTIMATE, IMPROVE>
      accum_i(estimates, lattice, estimator, fragments);
    for (int s = 0; s < nvs; ++s) {
      weight.start_bottom(s, time_t(0), s, spins[s]);
      accum_i.start_bottom(s, time_t(0), s, spins[s]);
      weight.stop_top(s, time_t(1), s, spins[s]);
      accum_i.stop_top(s, time_t(1), s, spins[s]);
    }
  }

  // accumulate cluster properties
  coll.set_num_clusters(nc);
  if (IMPROVE()) for (int c = 0; c < nc; ++c) coll += estimates[c];

  // determine whether clusters are flipped or not
  for (int c = 0; c < nc; ++c)
    estimates[c].to_flip = ((2*uniform_01()-1) < (FIELD() ? std::tanh(clusters[c].weight) : 0));

  // flip spins
  for (int s = 0; s < nvs; ++s) if (estimates[fragments[s].id()].to_flip) spins[s] ^= 1;

  //
  // measurement
  //

  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Volume"] << (double)lattice.volume();
  obs["Number of Sites"] << (double)nrs;
  obs["Number of Clusters"] << (double)coll.num_clusters();
  double ene = model.energy_offset() - nop / beta;
  if (model.has_field())
    for (int c = 0; c < nc; ++c)
      ene += (estimates[c].to_flip ? -clusters[c].weight : clusters[c].weight);
  coll.set_energy(ene);

  coll.commit(obs, lattice, beta, 1, nop);
}

typedef looper::evaluator<loop_config::measurement_set> loop_evaluator;

//
// dynamic registration to the factories
//

PARAPACK_REGISTER_ALGORITHM(loop_worker, "ising");
PARAPACK_REGISTER_ALGORITHM(alps::parapack::single_exchange_worker<loop_worker>, "ising exchange");

} // end namespace

#endif // LOOPER_ONLY_PATH_INTEGRAL
