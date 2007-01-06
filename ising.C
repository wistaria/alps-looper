/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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
#include "loop_factory.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/montecarlo.h>
#include <looper/permutation.h>
#include <looper/temperature.h>

namespace {

struct dummy_operator {
  static bool is_offdiagonal() { return false; }
  static int pos() { return 0; }
  static bool is_site() { return false; }
  static bool is_bond() { return false; }
};

class loop_worker {
public:
  typedef looper::classical mc_type;

  typedef looper::lattice_helper<loop_config::lattice_graph_t> lattice_t;
  typedef loop_config::loop_graph_t loop_graph_t;
  typedef loop_config::model_t model_t;
  typedef std::vector<dummy_operator> operator_string_t;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef looper::estimator<loop_config::measurement_set, mc_type, lattice_t, time_t>::type
    estimator_t;

  loop_worker(alps::Parameters const& p, alps::ObservableSet& obs);
  template<typename ENGINE>
  void run(ENGINE& eng, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void save(alps::ODump& dp) const { dp << mcs << spins; }
  void load(alps::IDump& dp) { dp >> mcs >> spins; }

protected:
  template<typename ENGINE>
  void build(ENGINE& eng);

  template<typename ENGINE, typename FIELD, typename IMPROVE>
  void flip(ENGINE& eng, alps::ObservableSet& obs);

  void measure(alps::ObservableSet& obs);

private:
  // helpers
  lattice_t lattice;
  model_t model;

  // parameters
  looper::temperature temperature;
  double beta;
  bool use_improved_estimator;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;

  estimator_t estimator;

  // working vectors
  int nop;
  std::vector<cluster_fragment_t> fragments;
  std::vector<cluster_info_t> clusters;
  std::vector<looper::estimate<estimator_t>::type> estimates;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::Parameters const& p, alps::ObservableSet& obs) :
  lattice(p), model(p, lattice), temperature(p), mcs(p) {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  use_improved_estimator = !model.has_field() && !p.defined("DISABLE_IMPROVED_ESTIMATOR");

  if (model.is_quantal())
    boost::throw_exception(std::invalid_argument("not classical Ising model"));
  if (model.has_field())
    std::cerr << "WARNING: model has magnetic field\n";
  if (model.is_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";
  if (!use_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  perm.resize(max_virtual_sites(lattice));

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);

  // measurements
  obs << make_observable(alps::SimpleRealObservable("Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Volume"));
  obs << make_observable(alps::SimpleRealObservable("Number of Sites"));
  obs << make_observable(alps::SimpleRealObservable("Number of Clusters"));
  looper::energy_estimator::initialize(obs, false);
  estimator.initialize(obs, p, lattice, false, use_improved_estimator);
}

template<typename ENGINE>
void loop_worker::run(ENGINE& eng, alps::ObservableSet& obs) {
  if (!mcs.can_work()) return;
  ++mcs;
  beta = 1.0 / temperature(mcs());

  build(eng);

  //           FIELD               IMPROVE
  flip<ENGINE, boost::mpl::true_,  boost::mpl::true_ >(eng, obs);
  flip<ENGINE, boost::mpl::true_,  boost::mpl::false_>(eng, obs);
  flip<ENGINE, boost::mpl::false_, boost::mpl::true_ >(eng, obs);
  flip<ENGINE, boost::mpl::false_, boost::mpl::false_>(eng, obs);

  measure(obs);
}


//
// cluster construction
//

template<typename ENGINE>
void loop_worker::build(ENGINE& eng) {
  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);

  // building up clusters
  nop = 0;
  boost::variate_generator<ENGINE&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());
  boost::variate_generator<ENGINE&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(beta * model.graph_weight()));
  double t = r_time();
  while (t < 1) {
    loop_graph_t g = model.choose_graph(r_uniform);
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
    BOOST_FOREACH(looper::real_site_descriptor<lattice_t>::type rs, sites(lattice.rg())) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      boost::tie(vsi, vsi_end) = sites(lattice, rs);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins.begin() + offset, r_uniform);
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, offset+perm[i]);
    }
  }
}


//
// cluster flip
//

template<typename ENGINE, typename FIELD, typename IMPROVE>
void loop_worker::flip(ENGINE& eng, alps::ObservableSet& obs) {
  if (FIELD() || use_improved_estimator != IMPROVE()) return;

  boost::variate_generator<ENGINE&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());

  int nvs = num_sites(lattice.vg());

  // assign cluster id
  int nc = 0;
  BOOST_FOREACH(cluster_fragment_t& f, fragments) if (f.is_root()) f.id = nc++;
  BOOST_FOREACH(cluster_fragment_t& f, fragments) f.id = cluster_id(fragments, f);

  clusters.resize(0); clusters.resize(nc);

  cluster_info_t::accumulator<cluster_fragment_t, FIELD,
    boost::mpl::false_, IMPROVE> weight(clusters, fragments, model.field(), 0, 0);
  looper::accumulator<estimator_t, cluster_fragment_t, IMPROVE>
    accum(estimates, nc, lattice, estimator, fragments);
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.at_bot(s, time_t(0), s, spins[s]);
    weight.at_top(s, time_t(1), s, spins[s]);
    accum.at_bot(s, time_t(0), s, spins[s]);
    accum.at_top(s, time_t(1), s, spins[s]);
  }

  // determine whether clusters are flipped or not
  BOOST_FOREACH(cluster_info_t& ci, clusters)
    ci.to_flip = ((2*r_uniform()-1) < (FIELD() ? std::tanh(ci.weight) : 0));

  // improved measurement
  if (IMPROVE()) {
    typename looper::collector<estimator_t>::type coll = get_collector(estimator);
    coll = std::accumulate(estimates.begin(), estimates.end(), coll);
    estimator.improved_measurement(obs, lattice, beta, 1, spins, operator_string_t(),
      spins, fragments, coll);
  }
  obs["Number of Clusters"] << (double)clusters.size();

  // flip spins
  for (int s = 0; s < nvs; ++s) if (clusters[fragments[s].id].to_flip) spins[s] ^= 1;
}


//
// measurement
//

void loop_worker::measure(alps::ObservableSet& obs) {
  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Volume"] << (double)lattice.volume();
  obs["Number of Sites"] << (double)num_sites(lattice.rg());

  // energy
  double ene = model.energy_offset() - nop / beta;
  if (model.has_field())
    BOOST_FOREACH(const cluster_info_t& c, clusters)
      ene += (c.to_flip ? -c.weight : c.weight);
  looper::energy_estimator::measurement(obs, lattice, beta, nop, 1, ene);

  estimator.normal_measurement(obs, lattice, beta, 1, spins, operator_string_t(), spins);
}


//
// dynamic registration to the factories
//

const bool loop_registered =
  loop_factory::instance()->register_worker<loop_worker>("Ising");
const bool evaluator_registered = evaluator_factory::instance()->
  register_evaluator<looper::evaluator<loop_config::measurement_set> >("Ising");

} // end namespace
