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
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/type.h>
#include <looper/union_find.h>

#ifndef LOOPER_ONLY_PATH_INTEGRAL

namespace {

class loop_worker {
public:
  typedef looper::sse mc_type;

  typedef looper::lattice_helper<loop_config::lattice_graph_t> lattice_t;
  typedef int time_t;
  typedef loop_config::loop_graph_t loop_graph_t;
  typedef loop_config::model_t model_t;
  typedef looper::local_operator<mc_type, loop_graph_t, time_t> local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef operator_string_t::iterator operator_iterator;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef looper::estimator<loop_config::measurement_set, mc_type, lattice_t, time_t>::type
    estimator_t;

  loop_worker(alps::Parameters const& p, alps::ObservableSet& obs);
  template<typename ENGINE>
  void run(ENGINE& eng, alps::ObservableSet& obs);

  void set_beta(double b) { temperature.set_beta(b); }
  double dlogw(double newbeta) const {
    return operators.size() * std::log(newbeta / beta);
  }

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void save(alps::ODump& dp) const { dp << mcs << spins << operators; }
  void load(alps::IDump& dp) { dp >> mcs >> spins >> operators; }

protected:
  template<typename ENGINE>
  void build(ENGINE& eng);

  template<typename ENGINE, typename FIELD, typename SIGN, typename IMPROVE>
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

loop_worker::loop_worker(alps::Parameters const& p, alps::ObservableSet& obs) :
  lattice(p), model(p, lattice, looper::is_path_integral<mc_type>::type()), temperature(p),
  mcs(p) {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  use_improved_estimator = !model.has_field() && !p.defined("DISABLE_IMPROVED_ESTIMATOR");

  if (model.has_field())
    boost::throw_exception(std::logic_error("longitudinal field is currently not supported "
      "in SSE representation"));
  if (model.is_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";
  if (model.is_signed())
    std::cerr << "WARNING: model has negative signs\n";
  if (!use_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  perm.resize(max_virtual_sites(lattice));

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);

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
  int nop = operators.size();
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  boost::variate_generator<ENGINE&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());
  double bw = beta * model.graph_weight();
  bool try_gap = true;
  for (operator_iterator opi = operators_p.begin(); try_gap || opi != operators_p.end();) {

    // diagonal update & labeling
    if (try_gap) {
      if ((nop+1) * r_uniform() < bw) {
        loop_graph_t g = model.choose_graph(r_uniform);
        if ((is_bond(g) && is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                                            spins_c[target(pos(g), lattice.vg())])) ||
            (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
          operators.push_back(local_operator_t(g));
          ++nop;
        } else {
          try_gap = false;
          continue;
        }
      } else {
        try_gap = false;
        continue;
      }
    } else {
      if (opi->is_diagonal()) {
        if (bw * r_uniform() < nop) {
          --nop;
          ++opi;
          continue;
        } else {
          if (opi->is_site()) {
            opi->assign_graph(model.choose_diagonal(r_uniform, opi->loc(), spins_c[opi->pos()]));
          } else {
            opi->assign_graph(model.choose_diagonal(r_uniform, opi->loc(),
              spins_c[source(opi->pos(), lattice.vg())],
              spins_c[target(opi->pos(), lattice.vg())]));
          }
        }
      } else {
        if (opi->is_bond())
          opi->assign_graph(model.choose_offdiagonal(r_uniform, opi->loc(),
            spins_c[source(opi->pos(), lattice.vg())],
            spins_c[target(opi->pos(), lattice.vg())]));
      }
      operators.push_back(*opi);
      ++opi;
      try_gap = true;
    }

    operator_iterator oi = operators.end() - 1;
    if (oi->is_bond()) {
      int s0 = source(oi->pos(), lattice.vg());
      int s1 = target(oi->pos(), lattice.vg());
      if (oi->is_offdiagonal()) {
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

  // symmetrize spins
  if (max_virtual_sites(lattice) == 1) {
    for (int i = 0; i < nvs; ++i) unify(fragments, i, current[i]);
  } else {
    BOOST_FOREACH(looper::real_site_descriptor<lattice_t>::type rs, sites(lattice.rg())) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      boost::tie(vsi, vsi_end) = sites(lattice, rs);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins_c.begin() + offset, r_uniform);
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, current[offset+perm[i]]);
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
  int nop = operators.size();

  // assign cluster id
  int nc = 0;
  BOOST_FOREACH(cluster_fragment_t& f, fragments) if (f.is_root()) f.id = nc++;
  BOOST_FOREACH(cluster_fragment_t& f, fragments) f.id = cluster_id(fragments, f);
  clusters.resize(0); clusters.resize(nc);

  std::copy(spins.begin(), spins.end(), spins_c.begin());
  cluster_info_t::accumulator<cluster_fragment_t, FIELD, SIGN, IMPROVE>
    weight(clusters, fragments, model.field(), model.bond_sign(), model.site_sign());
  looper::accumulator<estimator_t, cluster_fragment_t, IMPROVE>
    accum(estimates, nc, lattice, estimator, fragments);
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.at_bot(s, time_t(0), s, spins_c[s]);
    accum.at_bot(s, time_t(0), s, spins_c[s]);
  }
  int t = 0;
  int negop = 0; // number of operators with negative weights
  BOOST_FOREACH(local_operator_t& op, operators) {
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
    ++t;
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.at_top(current[s], time_t(nop), s, spins_c[s]);
    accum.at_top(current[s], time_t(nop), s, spins_c[s]);
  }
  sign = ((negop & 1) == 1) ? -1 : 1;

  // determine whether clusters are flipped or not
  double improved_sign = sign;
  BOOST_FOREACH(cluster_info_t& ci, clusters) {
    ci.to_flip = ((2*r_uniform()-1) < (FIELD() ? std::tanh(beta * ci.weight) : 0));
    if (SIGN() && IMPROVE() && (ci.sign & 1 == 1)) improved_sign = 0;
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
    if (clusters[fragments[op.loop_0()].id].to_flip ^ clusters[fragments[op.loop_1()].id].to_flip)
      op.flip();
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id].to_flip) spins[s] ^= 1;
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
  looper::energy_estimator::measurement(obs, lattice, beta, nop, sign, ene);

  // other quantities
  estimator.normal_measurement(obs, lattice, beta, sign, spins, operators, spins_c);
}


//
// dynamic registration to the factories
//

const bool worker_registered =
  loop_factory::instance()->register_worker<loop_worker>("SSE");
const bool evaluator_registered = loop_factory::instance()->
  register_evaluator<looper::evaluator<loop_config::measurement_set> >("SSE");

} // end namespace

#endif // LOOPER_ONLY_PATH_INTEGRAL
