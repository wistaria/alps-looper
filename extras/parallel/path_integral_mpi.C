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

// #define PARA_RANGE_DEBUG 1

#include "loop_config.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/type.h>
#include <looper/union_find.h>
#include <looper/parallel.h>
#include <alps/parapack/parallel.h>

namespace {

namespace mpi = boost::mpi;

class loop_worker : public alps::parapack::mc_worker, private loop_config {
public:
  typedef looper::path_integral mc_type;

  typedef looper::local_operator<mc_type, loop_graph_t, time_t> local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef operator_string_t::iterator operator_iterator;

  typedef looper::union_find::node cluster_fragment_t;
  typedef looper::cluster_info cluster_info_t;

  typedef looper::estimator<measurement_set, mc_type, lattice_t, time_t>::type estimator_t;

  typedef looper::parallel_cluster_unifier<looper::estimate<estimator_t>::type,
    looper::collector<estimator_t>::type> unifier_t;

  loop_worker(mpi::communicator const& c, alps::Parameters const& p);
  virtual ~loop_worker() {}

  void init_observables(alps::Parameters const& params, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs);

  void save(alps::ODump& dp) const { dp << mcs << spins << spins_t << operators; }
  void load(alps::IDump& dp) { dp >> mcs >> spins >> spins_t >> operators; }

protected:
  void build();

  template<typename FIELD, typename SIGN, typename IMPROVE>
  void flip(alps::ObservableSet& obs);

  boost::tuple<int, int> real_site_range(int nrs, int nprocs, int rank) const;

private:
  // helpers
  mpi::communicator comm;
  lattice_t lattice;
  model_t model;
  unifier_t unifier;

  // parameters
  looper::temperature temperature;
  double beta;
  double tau0, tau1;
  bool use_improved_estimator;

  // index for parallel
  int rs_local_begin, rs_local_end;
  int vs_local_begin, vs_local_end, nvs_local;

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
  std::vector<unifier_t::flip_t> to_flip;
  std::vector<looper::estimate<estimator_t>::type> estimates;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(mpi::communicator const& c, alps::Parameters const& p)
  : alps::parapack::mc_worker(p), comm(c), lattice(p),
    model(p, lattice, /* is_path_integral = */ true),
    unifier(comm, num_sites(lattice.vg())),
    temperature(p), mcs(p) {
  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(/* support_longitudinal_field = */ false,
    /* support_negative_sign = */ false);
  use_improved_estimator = (!model.has_field()) && (!p.defined("DISABLE_IMPROVED_ESTIMATOR"));
  if (!use_improved_estimator) std::cerr << "WARNING: improved estimator is disabled\n";

  tau0 = 1.0 * comm.rank() / comm.size();
  tau1 = 1.0 * (comm.rank() + 1) / comm.size();

  // configuration
  int nrs = num_sites(lattice.rg());
  int nvs = num_sites(lattice.vg());
  boost::tie(rs_local_begin, rs_local_end) = real_site_range(nrs, comm.size(), comm.rank());
  looper::virtual_site_iterator<lattice_t>::type vs, ve;
  boost::tie(vs, ve) = sites(lattice, rs_local_begin); vs_local_begin = *vs;
  boost::tie(vs, ve) = sites(lattice, rs_local_end-1); vs_local_end = *ve;
  nvs_local = vs_local_end - vs_local_begin;
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_t.resize(nvs_local); std::fill(spins_t.begin(), spins_t.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);
  perm.resize(max_virtual_sites(lattice));

  // initialize estimators
  estimator.initialize(p, lattice, model.is_signed(), use_improved_estimator);
}

void loop_worker::init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
  if (comm.rank() == 0) {
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
  looper::energy_estimator::init_observables(obs, model.is_signed());
  estimator.init_observables(obs, model.is_signed());
}

void loop_worker::run(alps::ObservableSet& obs) {
  // if (!mcs.can_work()) return;
  ++mcs;

  beta = 1.0 / temperature(mcs());
  tau0 = 1.0 * comm.rank() / comm.size();
  tau1 = 1.0 * (comm.rank() + 1) / comm.size();

  build();

  //   FIELD               SIGN                IMPROVE
  flip<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs);
  flip<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs);
}


//
// diagonal update and cluster construction
//

void loop_worker::build() {
  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(2 * nvs + nvs_local);
  for (int s = 0; s < nvs; ++s) current[s] = s;
  int offset = 2 * nvs - vs_local_begin;
  for (int s = vs_local_begin; s < vs_local_end; ++s) current[s] = offset + s;

  boost::variate_generator<engine_type&, boost::exponential_distribution<> >
    r_time(engine(), boost::exponential_distribution<>(beta * model.graph_weight()));
  double t = tau0 + r_time();
  for (operator_iterator opi = operators_p.begin(); t < tau1 || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      loop_graph_t g = model.choose_graph(generator_01());
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
        oi->assign_graph(model.choose_offdiagonal(generator_01(), oi->loc(),
          spins_c[s0], spins_c[s1]));
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
  if (max_virtual_sites(lattice) == 1) {
    for (int s = vs_local_begin; s < vs_local_end; ++s) unify(fragments, s, 2*nvs - vs_local_begin + s);
  } else {
    for (int rs = rs_local_begin; rs < rs_local_end; ++rs) {
      looper::virtual_site_iterator<lattice_t>::type vsi_sta, vsi_end;
      boost::tie(vsi_sta, vsi_end) = sites(lattice, rs);
      int offset1 = *vsi_sta;
      int offset2 = *vsi_sta - vs_local_begin;
      int s2 = *vsi_end - *vsi_sta;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins_t.begin() + offset2, spins.begin() + offset1, generator_01());
      for (int i = 0; i < s2; ++i) unify(fragments, offset1+i, 2*nvs+offset2+perm[i]);
    }
  }

  for (int s = 2 * nvs - 1; s >= 0; --s) set_root(fragments, s);
}


//
// cluster flip
//

template<typename FIELD, typename SIGN, typename IMPROVE>
void loop_worker::flip(alps::ObservableSet& obs) {
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      use_improved_estimator != IMPROVE()) return;

  int nvs = num_sites(lattice.vg());

  // assign cluster id
  int nc = 0;
  for (int s = 0; s < 2 * nvs; ++s) if (fragments[s].is_root()) fragments[s].set_id(nc++);
  int noc = nc;
  for (int s = 2 * nvs; s < fragments.size(); ++s)
    if (fragments[s].is_root()) fragments[s].set_id(nc++);
  BOOST_FOREACH(cluster_fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
  to_flip.resize(nc);

  std::copy(spins.begin(), spins.end(), spins_c.begin());
  looper::accumulator<estimator_t, cluster_fragment_t, IMPROVE>
    accum(estimates, nc, lattice, estimator, fragments);
  if (comm.rank() == 0) {
    for (unsigned int s = 0; s < vs_local_begin; ++s)
      accum.start_bottom(s, time_t(tau0), s, spins_c[s]);
    for (unsigned int s = vs_local_begin; s < vs_local_end; ++s)
      accum.start_bottom(2*nvs+s-vs_local_begin, time_t(tau0), s, spins_c[s]);
    for (unsigned int s = vs_local_end; s < nvs; ++s)
      accum.start_bottom(s, time_t(tau0), s, spins_c[s]);
  } else {
    for (unsigned int s = 0; s < vs_local_begin; ++s)
      accum.start(s, time_t(tau0), s, spins_c[s]);
    for (unsigned int s = vs_local_begin; s < vs_local_end; ++s)
      accum.start(2*nvs+s-vs_local_begin, time_t(tau0), s, spins_c[s]);
    for (unsigned int s = vs_local_end; s < nvs; ++s)
      accum.start(s, time_t(tau0), s, spins_c[s]);
  }
  BOOST_FOREACH(local_operator_t& op, operators) {
    time_t t = op.time();
    if (op.is_bond()) {
      if (!op.is_frozen_bond_graph()) {
        int b = op.pos();
        int s0 = source(b, lattice.vg());
        int s1 = target(b, lattice.vg());
        accum.end_b(op.loop_l0(), op.loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
        if (op.is_offdiagonal()) {
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
        }
        accum.begin_b(op.loop_u0(), op.loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
      }
    } else {
      if (!op.is_frozen_site_graph()) {
        int s = op.pos();
        accum.end_s(op.loop_l(), t, s, spins_c[s]);
        if (op.is_offdiagonal()) {
          spins_c[s] ^= 1;
        }
        accum.begin_s(op.loop_u(), t, s, spins_c[s]);
      }
    }
  }
  if (comm.rank() == (comm.size() - 1)) {
    for (unsigned int s = 0; s < nvs; ++s) accum.stop_top(current[s], time_t(tau1), s, spins_c[s]);
  } else {
    for (unsigned int s = 0; s < nvs; ++s) accum.stop(current[s], time_t(tau1), s, spins_c[s]);
  }

  // accumulate cluster properties
  typename looper::collector<estimator_t>::type coll = get_collector(estimator);
  coll.set_num_operators(operators.size());
  coll.set_num_open_clusters(noc);
  coll.set_num_clusters(nc - noc);
  for (unsigned int c = noc; c < nc; ++c) coll += estimates[c];

  // determine whether clusters are flipped or not
  unifier.unify(coll, to_flip, fragments, estimates, generator_01());
  for (int c = noc; c < nc; ++c) to_flip[c].set_flip(uniform_01() < 0.5);

  // improved measurement
  if (comm.rank() == 0)
    estimator.improved_measurement(obs, lattice, beta, 1, spins, operators, spins_c,
      fragments, coll);

  // flip operators & spins
  BOOST_FOREACH(local_operator_t& op, operators)
    if (to_flip[fragments[op.loop_0()].id()] ^ to_flip[fragments[op.loop_1()].id()]) op.flip();
  for (int s = 0; s < vs_local_begin; ++s) if (to_flip[fragments[s].id()]) spins[s] ^= 1;
  int offset1 = vs_local_begin;
  for (int s = vs_local_begin; s < vs_local_end; ++s) {
    if (to_flip[fragments[s].id()]) spins_t[s - offset1] ^= 1;
    if (to_flip[fragments[2*nvs + s - offset1].id()]) spins[s] ^= 1;
  }
  for (int s = vs_local_end; s < nvs; ++s) if (to_flip[fragments[s].id()]) spins[s] ^= 1;

  // measurement
  if (comm.rank() == 0) {
    obs["Temperature"] << 1/beta;
    obs["Inverse Temperature"] << beta;
    obs["Volume"] << (double)lattice.volume();
    obs["Number of Sites"] << (double)num_sites(lattice.rg());
    obs["Number of Clusters"] << coll.num_clusters();

    double nop = coll.num_operators();
    double ene = model.energy_offset() - nop / beta;
    looper::energy_estimator::measurement(obs, lattice, beta, nop, 1, ene);
    estimator.normal_measurement(obs, lattice, beta, 1, spins, operators, spins_c);
  }
}

boost::tuple<int, int> loop_worker::real_site_range(int nrs, int nprocs, int rank) const {
  int nrs_local = nrs / nprocs;
  int nextra = nrs % nprocs;
  int rs_begin = rank * nrs_local + (rank < nextra ? rank : nextra);
  int rs_local_end = rs_begin + nrs_local + (rank < nextra ? 1 : 0);
#ifdef PARA_RANGE_DEBUG
  comm.barrier();
  if (rank == 0)
    std::cerr << "rank " << rank << ": real_site_range: nrs = " << nrs << std::endl;
  comm.barrier();
  std::cerr << nrs_local << ' ' << nextra << std::endl;
  std::cerr << "rank " << rank << ": real_site_range: rs_begin = " << rs_begin
            << ", rs_local_end = " << rs_local_end << std::endl;
  comm.barrier();
#endif
  return boost::make_tuple(rs_begin, rs_local_end);
}

typedef looper::evaluator<loop_config::measurement_set> loop_evaluator;

//
// dynamic registration to the factories
//

PARAPACK_REGISTER_PARALLEL_WORKER(loop_worker, "path integral");
PARAPACK_REGISTER_EVALUATOR(loop_evaluator, "path integral");

} // end namespace
