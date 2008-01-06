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

#include "loop_config.h"
#include "loop_factory.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/type.h>
#include <looper/union_find.h>
#include <looper/parallel.h>

namespace {

class loop_worker {
public:
  typedef looper::path_integral mc_type;

  typedef looper::lattice_helper<loop_config::lattice_graph_t> lattice_t;
  typedef loop_config::time_t time_t;
  typedef loop_config::loop_graph_t loop_graph_t;
  typedef loop_config::model_t model_t;
  typedef looper::local_operator<mc_type, loop_graph_t, time_t> local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef operator_string_t::iterator operator_iterator;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef looper::estimator<loop_config::measurement_set, mc_type, lattice_t, time_t>::type
    estimator_t;

  loop_worker(alps::communicator_helper const& c, alps::Parameters const& p,
    alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  template<typename ENGINE>
  void run(ENGINE& eng, alps::ObservableSet& obs);

  // for exchange Monte Carlo
  void set_beta(double b) { temperature.set_beta(b); }
  double g_weight() const { return operators.size(); }
  double lambda(double beta) const { return std::log(beta); }

  void save(alps::ODump& dp) const { dp << mcs << spins << spins_t << operators; }
  void load(alps::IDump& dp) { dp >> mcs >> spins >> spins_t >> operators; }

protected:
  template<typename ENGINE>
  void build(ENGINE& eng);

  template<typename ENGINE, typename FIELD, typename SIGN, typename IMPROVE>
  void flip(ENGINE& eng, alps::ObservableSet& obs);

  void measure(alps::ObservableSet& obs);

private:
  // helpers
  alps::communicator_helper comm;
  lattice_t lattice;
  model_t model;
  looper::parallel_cluster_unifier unifier;

  // parameters
  looper::temperature temperature;
  double beta;
  double tau0, tau1;
  bool use_improved_estimator;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<int> spins_t; // only for master process
  std::vector<local_operator_t> operators;

  // observables
  double sign;
  estimator_t estimator;

  // working vectors
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<cluster_info_t> clusters;
  std::vector<looper::estimate<estimator_t>::type> estimates;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::communicator_helper const& c, alps::Parameters const& p,
  alps::ObservableSet& obs) :
  comm(c), lattice(p), model(p, lattice, looper::is_path_integral<mc_type>::type()),
  unifier(c, num_sites(lattice.vg())), temperature(p), mcs(p) {
  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  use_improved_estimator = !model.has_field() && !p.defined("DISABLE_IMPROVED_ESTIMATOR");

  if (model.has_field())
    std::cerr << "WARNING: model has a magnetic field\n";
  if (model.is_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";
  if (model.is_signed())
    std::cerr << "WARNING: model has negative signs\n";
  if (!use_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  if (is_master(comm))
    perm.resize(max_virtual_sites(lattice));

  tau0 = 1.0 * process_id(comm) / num_processes(comm);
  tau1 = 1.0 * (process_id(comm) + 1) / num_processes(comm);

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  if (is_master(comm)) spins_t.resize(nvs); std::copy(spins.begin(), spins.end(), spins_t.begin());
  spins_c.resize(nvs);
  current.resize(nvs);

  if (is_master(comm)) {
    obs << make_observable(alps::SimpleRealObservable("Temperature"));
    obs << make_observable(alps::SimpleRealObservable("Inverse Temperature"));
    obs << make_observable(alps::SimpleRealObservable("Volume"));
    obs << make_observable(alps::SimpleRealObservable("Number of Sites"));
    obs << make_observable(alps::SimpleRealObservable("Number of Clusters"));
    if (model.is_signed()) {
      obs << alps::RealObservable("Sign");
      if (use_improved_estimator) {
        obs << alps::RealObservable("Weight of Zero-Meron Sector");
        obs << alps::RealObservable("Sign in Zero-Meron Sector");
      }
    }
  }
  looper::energy_estimator::initialize(obs, model.is_signed());
  estimator.initialize(obs, p, lattice, model.is_signed(), use_improved_estimator);
}

template<typename ENGINE>
void loop_worker::run(ENGINE& eng, alps::ObservableSet& obs) {
  // if (!mcs.can_work()) return;
  ++mcs;
  beta = 1.0 / temperature(mcs());

  build(eng);

  //           FIELD               SIGN                IMPROVE
  flip<ENGINE, boost::mpl::true_,  boost::mpl::true_,  boost::mpl::true_ >(eng, obs);
  flip<ENGINE, boost::mpl::true_,  boost::mpl::true_,  boost::mpl::false_>(eng, obs);
  flip<ENGINE, boost::mpl::true_,  boost::mpl::false_, boost::mpl::true_ >(eng, obs);
  flip<ENGINE, boost::mpl::true_,  boost::mpl::false_, boost::mpl::false_>(eng, obs);
  flip<ENGINE, boost::mpl::false_, boost::mpl::true_,  boost::mpl::true_ >(eng, obs);
  flip<ENGINE, boost::mpl::false_, boost::mpl::true_,  boost::mpl::false_>(eng, obs);
  flip<ENGINE, boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(eng, obs);
  flip<ENGINE, boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(eng, obs);

  measure(obs);
}


//
// diagonal update and cluster construction
//

template<typename ENGINE>
void loop_worker::build(ENGINE& eng) {
  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0);
  if (is_master(comm)) {
    fragments.resize(3 * nvs);
    for (int s = 0; s < nvs; ++s) current[s] = 2 * nvs + s;
  } else {
    fragments.resize(2 * nvs);
    for (int s = 0; s < nvs; ++s) current[s] = s;
  }

  boost::variate_generator<ENGINE&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());
  boost::variate_generator<ENGINE&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(beta * model.graph_weight()));
  double t = tau0 + r_time();
  for (operator_iterator opi = operators_p.begin(); t < tau1 || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      loop_graph_t g = model.choose_graph(r_uniform);
      if ((is_bond(g) && is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                                          spins_c[target(pos(g), lattice.vg())])) ||
          (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
        operators.push_back(local_operator_t(g, t));
        t += r_time();
      } else {
        t += r_time();
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
      int s0 = source(oi->pos(), lattice.vg());
      int s1 = target(oi->pos(), lattice.vg());
      if (oi->is_offdiagonal()) {
        oi->assign_graph(model.choose_offdiagonal(r_uniform, oi->loc(), spins_c[s0], spins_c[s1]));
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
      boost::tie(current[s0], current[s1], oi->loop0, oi->loop1) =
        reconnect(fragments, oi->graph(), current[s0], current[s1]);
    } else {
      int s = oi->pos();
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
      boost::tie(current[s], oi->loop0, oi->loop1) = reconnect(fragments, oi->graph(), current[s]);
    }
  }

  // connect to top
  for (int s = 0; s < nvs; ++s) unify(fragments, current[s], nvs + s);

  // symmetrize spins
  if (is_master(comm)) {
    if (max_virtual_sites(lattice) == 1) {
      for (int s = 0; s < nvs; ++s) unify(fragments, s, 2*nvs + s);
    } else {
      BOOST_FOREACH(looper::real_site_descriptor<lattice_t>::type rs, sites(lattice.rg())) {
        looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
        boost::tie(vsi, vsi_end) = sites(lattice, rs);
        int offset = *vsi;
        int s2 = *vsi_end - *vsi;
        for (int i = 0; i < s2; ++i) perm[i] = i;
        looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
          spins_t.begin() + offset, spins.begin() + offset, r_uniform);
        for (int i = 0; i < s2; ++i) unify(fragments, offset+i, 2*nvs+offset+perm[i]);
      }
    }
  }
}


//
// cluster flip
//

template<typename ENGINE, typename FIELD, typename SIGN, typename IMPROVE>
void loop_worker::flip(ENGINE& eng, alps::ObservableSet& obs) {
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      use_improved_estimator != IMPROVE()) return;

  boost::variate_generator<ENGINE&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());

  int nvs = num_sites(lattice.vg());

  // assign cluster id
  for (int v = 0; v < 2*nvs; ++v) {
    int r = root_index(fragments, v);
    if (r > v) {
      fragments[v].set_as_root(fragments[r].weight());
      fragments[r].set_parent(v);
    }
  }
  int nc = 0;
  for (int v = 0; v < 2*nvs; ++v) if (fragments[v].is_root()) fragments[v].set_id(nc++);
  int noc = nc;
  for (int v = 2*nvs; v < fragments.size(); ++v)
    if (fragments[v].is_root()) fragments[v].set_id(nc++);
  BOOST_FOREACH(cluster_fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
  clusters.resize(0); clusters.resize(nc);

  std::copy(spins.begin(), spins.end(), spins_c.begin());
  cluster_info_t::accumulator<cluster_fragment_t, FIELD, SIGN, IMPROVE>
    weight(clusters, fragments, model.field(), model.bond_sign(), model.site_sign());
  looper::accumulator<estimator_t, cluster_fragment_t, IMPROVE>
    accum(estimates, nc, lattice, estimator, fragments);
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.at_bot(s, time_t(tau0), s, spins_c[s]);
    accum.at_bot(s, time_t(tau0), s, spins_c[s]);
  }
  int negop = 0; // number of operators with negative weights
  BOOST_FOREACH(local_operator_t& op, operators) {
    time_t t = op.time();
    if (op.is_bond()) {
      if (!op.is_frozen_bond_graph()) {
        int b = op.pos();
        int s0 = source(b, lattice.vg());
        int s1 = target(b, lattice.vg());
        weight.term_b(op.loop_l0(), op.loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        accum.term_b(op.loop_l0(), op.loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        if (op.is_offdiagonal()) {
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
          if (SIGN()) negop += model.bond_sign(op.pos());
        }
        weight.start_b(op.loop_u0(), op.loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        accum.start_b(op.loop_u0(), op.loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
      }
    } else {
      if (!op.is_frozen_site_graph()) {
        int s = op.pos();
        weight.term_s(op.loop_l(), t, s, spins_c[s]);
        accum.term_s(op.loop_l(), t, s, spins_c[s]);
        if (op.is_offdiagonal()) {
          spins_c[s] ^= 1;
          if (SIGN()) negop += model.site_sign(op.pos());
        }
        weight.start_s(op.loop_u(), t, s, spins_c[s]);
        accum.start_s(op.loop_u(), t, s, spins_c[s]);
      }
    }
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.at_top(current[s], time_t(tau1), s, spins_c[s]);
    accum.at_top(current[s], time_t(tau1), s, spins_c[s]);
  }
  sign = ((negop & 1) == 1) ? -1 : 1;

  // global unification of clusters
  std::vector<looper::parallel_cluster_unifier::flip_t> const& to_flip =
    unifier.unify(noc, fragments, r_uniform);
  for (int c = 0; c < noc; ++c) clusters[c].to_flip = to_flip[c].flip();

  // determine whether clusters are flipped or not
  double improved_sign = sign;
  for (int c = noc; c < nc; ++c) {
    clusters[c].to_flip =
      ((2*r_uniform()-1) < (FIELD() ? std::tanh(beta * clusters[c].weight) : 0));
    if (SIGN() && IMPROVE() && (clusters[c].sign & 1 == 1)) improved_sign = 0;
  }

  // improved measurement
  if (IMPROVE()) {
    typename looper::collector<estimator_t>::type coll = get_collector(estimator);
    coll = std::accumulate(estimates.begin(), estimates.end(), coll);
    estimator.improved_measurement(obs, lattice, beta, improved_sign, spins, operators,
      spins_c, fragments, coll);
    if (SIGN()) {
      obs["Sign"] << improved_sign;
      if (alps::is_zero(improved_sign)) {
        obs["Weight of Zero-Meron Sector"] << 0.;
      } else {
        obs["Weight of Zero-Meron Sector"] << 1.;
        obs["Sign in Zero-Meron Sector"] << improved_sign;
      }
    }
  }
  obs["Number of Clusters"] << (double)clusters.size();

  // flip operators & spins
  BOOST_FOREACH(local_operator_t& op, operators)
    if (clusters[fragments[op.loop_0()].id()].to_flip ^
        clusters[fragments[op.loop_1()].id()].to_flip) op.flip();
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id()].to_flip) spins[s] ^= 1;
  if (is_master(comm))
    for (int s = 0; s < nvs; ++s)
      if (clusters[fragments[2*nvs + s].id()].to_flip) spins_t[s] ^= 1;
}


//
// measurement
//

void loop_worker::measure(alps::ObservableSet& obs) {
  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Volume"] << (double)lattice.volume();
  obs["Number of Sites"] << (double)num_sites(lattice.rg());

  // sign
  if (model.is_signed() && !use_improved_estimator) obs["Sign"] << sign;

  // energy
  int nop = operators.size();
  double ene = model.energy_offset() - nop / beta;
  if (model.has_field())
    BOOST_FOREACH(const cluster_info_t& c, clusters)
      ene += (c.to_flip ? -c.weight : c.weight);
  looper::energy_estimator::measurement(obs, lattice, beta, nop, sign, ene);

  // other quantities
  estimator.normal_measurement(obs, lattice, beta, sign, spins, operators, spins_c);
}


//
// dynamic registration to the factories
//

// const bool worker_registered =
//   loop_factory::instance()->register_worker<loop_worker>("path integral");
// const bool evaluator_registered = loop_factory::instance()->
//   register_evaluator<looper::evaluator<loop_config::measurement_set> >("path integral");

} // end namespace