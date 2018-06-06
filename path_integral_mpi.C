/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>,
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

#include <alps/config.h>
#if defined(LOOPER_ENABLE_OPENMP) && defined(ALPS_ENABLE_OPENMP_WORKER) && !defined(LOOPER_OPENMP)
# define LOOPER_OPENMP
#endif

// #define ALPS_ENABLE_TIMER
// #define ALPS_TRACE_TIMER
// #define STD_OUTPUT 1

#include "loop_config.h"
#include <looper/capacity_mpi.h>
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/expand.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/timer_mpi.hpp>
#include <looper/type.h>
#include <looper/union_find.h>
#include <looper/parallel.h>
#include <alps/parapack/worker.h>

#ifdef LOOPER_OPENMP
# include <looper/padded_vector.h>
# include <looper/subaccumulate.h>
# include <omp.h>
#endif

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

  typedef boost::exponential_distribution<> expdist_t;

  loop_worker(mpi::communicator const& c, alps::Parameters const& p);
  virtual ~loop_worker() {
#ifdef LOOPER_OPENMP
    if (reserved) {
      looper::vector_capacity capacity(times_g, operators_g, operators_pg, estimates_ig, estimates_ng, fragments);
      capacity.report(comm);
    }
#else
    if (reserved) {
      looper::vector_capacity capacity(times, operators, operators_p, estimates_i, estimates_n, fragments);
      capacity.report(comm);
    }
#endif
    if (enable_improved_estimator) {
      unifier_i.capacity_report();
    } else {
      unifier_n.capacity_report();
    }
    timer.stop(1);
    timer.summarize();
  }

  void init_observables(alps::Parameters const& params, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs);

  void save(alps::ODump& dp) const {
#ifdef LOOPER_OPENMP
    dp << mcs << spins << spins_t << operators_g;
#else
    dp << mcs << spins << spins_t << operators;
#endif
  }
  void load(alps::IDump& dp) {
#ifdef LOOPER_OPENMP
    dp >> mcs >> spins >> spins_t >> operators_g;
#else
    dp >> mcs >> spins >> spins_t >> operators;
#endif
  }

protected:
#ifdef LOOPER_OPENMP
  template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE,
    typename UNIFIER>
  void dispatch(alps::ObservableSet& obs, std::vector<COLLECTOR>& coll,
    std::vector<std::vector<ESTIMATE> >& estimates, UNIFIER& unifier);
#else
  template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE,
    typename UNIFIER>
  void dispatch(alps::ObservableSet& obs, COLLECTOR& coll, std::vector<ESTIMATE>& estimates,
    UNIFIER& unifier);
#endif

  boost::tuple<int, int> real_site_range(int nrs, int nprocs, int rank) const;

private:
  // helpers
  mpi::communicator comm;
  lattice_t lattice;
#ifdef LOOPER_OPENMP
  looper::simple_lattice_sharing sharing;
#endif
  model_t model;
  looper::parallel_cluster_unifier<estimator_t::improved_estimator::estimate,
    estimator_t::improved_estimator::collector> unifier_i;
  looper::parallel_cluster_unifier<estimator_t::normal_estimator::estimate,
    estimator_t::normal_estimator::collector> unifier_n;

  // parameters
  looper::temperature temperature;
  double beta;
  double tau0, tau1;
  bool enable_improved_estimator;

  // index for parallel
  int rs_local_begin, rs_local_end;
  int vs_local_begin, vs_local_end, nvs_local;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<int> spins_t;
#ifdef LOOPER_OPENMP
  std::vector<std::vector<local_operator_t> > operators_g;
#else
  std::vector<local_operator_t> operators;
#endif

  // observables
  estimator_t estimator;
#ifdef LOOPER_OPENMP
  std::vector<estimator_t::improved_estimator::collector> coll_ig;
  std::vector<estimator_t::normal_estimator::collector> coll_ng;
#else
  estimator_t::improved_estimator::collector coll_i;
  estimator_t::normal_estimator::collector coll_n;
#endif

  // working vectors
  std::vector<int> spins_c;
  std::vector<int> current;
  std::vector<cluster_fragment_t> fragments;
#ifdef LOOPER_OPENMP
  std::vector<std::vector<local_operator_t> > operators_pg;
  std::vector<std::vector<double> > times_g;
  std::vector<int> fragment_offset_g;
  std::vector<int> num_fragments_g;
  std::vector<int> noc_g;
  std::vector<int> ncc_g;
  std::vector<std::vector<estimator_t::improved_estimator::estimate> > estimates_ig;
  std::vector<std::vector<estimator_t::normal_estimator::estimate> > estimates_ng;
  std::vector<std::vector<int> > perm_g;
#else
  std::vector<local_operator_t> operators_p;
  std::vector<double> times;
  int fragment_offset;
  int num_fragments;
  std::vector<estimator_t::improved_estimator::estimate> estimates_i;
  std::vector<estimator_t::normal_estimator::estimate> estimates_n;
  std::vector<int> perm;
#endif

  alps::parapack::timer_mpi timer;
  bool reserved;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(mpi::communicator const& c, alps::Parameters const& p)
  : alps::parapack::mc_worker(p), comm(c), lattice(p),
#ifdef LOOPER_OPENMP
    sharing(/* will be initialized in model */),
    model(p, lattice, sharing, /* is_path_integral = */ true),
#else
    model(p, lattice, /* is_path_integral = */ true),
#endif
    unifier_i(comm), unifier_n(comm),
    temperature(p), mcs(p), timer(comm) {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(/* support_longitudinal_field = */ false,
                        /* support_negative_sign = */ false);

  enable_improved_estimator = (!model.has_field()) && (!p.defined("DISABLE_IMPROVED_ESTIMATOR"));
  if (!enable_improved_estimator && comm.rank() == 0)
    std::cout << "WARNING: improved estimator is disabled\n";

  tau0 = 1.0 * comm.rank() / comm.size();
  tau1 = 1.0 * (comm.rank() + 1) / comm.size();

  // configuration
#ifdef LOOPER_OPENMP
  int num_threads = omp_get_max_threads();
  if (comm.rank() == 0)
    std::cout << "Info: hybrid: " << c.size() << " x " << num_threads << std::endl;
#endif
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
#ifdef LOOPER_OPENMP
  operators_g.resize(num_threads);
  operators_pg.resize(num_threads);
  perm_g.resize(num_threads);
  for (int i = 0; i < num_threads; ++i) perm_g[i].resize(max_virtual_sites(lattice));
  coll_ig.resize(num_threads);
  coll_ng.resize(num_threads);
#else
  perm.resize(max_virtual_sites(lattice));
#endif

  // working vectors
  int reserve_times = p.value_or_default("RESERVE_TIMES", 0);
  int reserve_operators = p.value_or_default("RESERVE_OPERATORS", 0);
  int reserve_estimates = p.value_or_default("RESERVE_ESTIMATES", 0);
  int reserve_fragments = p.value_or_default("RESERVE_FRAGMENTS", 0);
  reserved = reserve_times || reserve_operators || reserve_estimates || reserve_fragments;
#ifdef LOOPER_OPENMP
  times_g.resize(num_threads);
  fragment_offset_g.resize(num_threads + 1);
  num_fragments_g.resize(num_threads);
  noc_g.resize(num_threads);
  ncc_g.resize(num_threads);
  estimates_ig.resize(num_threads);
  estimates_ng.resize(num_threads);
  if (reserved) {
    for (int tid = 0; tid < num_threads; ++tid) {
      times_g[tid].reserve(reserve_times);
      operators_g[tid].reserve(reserve_operators);
      operators_pg[tid].reserve(reserve_operators);
      if (enable_improved_estimator) {
        estimates_ig[tid].reserve(reserve_estimates);
      } else {
        estimates_ng[tid].reserve(reserve_estimates);
      }
    }
    fragments.reserve(reserve_fragments);
    looper::vector_capacity capacity(times_g, operators_g, operators_pg, estimates_ig,
      estimates_ng, fragments);
    capacity.report(comm);
  }
#else
  times.reserve(reserve_times);
  operators.reserve(reserve_operators);
  operators_p.reserve(reserve_operators);
  if (enable_improved_estimator) {
    estimates_i.reserve(reserve_estimates);
  } else {
    estimates_n.reserve(reserve_estimates);
  }
  fragments.reserve(reserve_fragments);
  looper::vector_capacity capacity(times, operators, operators_p, estimates_i, estimates_n,
    fragments);
  capacity.report(comm);
#endif

  // initialize timer
  timer.registrate( 1, "alps::parapack::scheduler::start");
  timer.registrate( 2, " loop_worker::run,all");
  timer.registrate( 3, "  dispatch,all");
  timer.registrate( 4, "   dispatch,init_spin&operator_info");
  timer.registrate( 5, "   dispatch,fill_times");
  timer.registrate( 6, "   dispatch,initialize_cluster_info");
  timer.registrate( 7, "   dispatch,insert&remove_operators");
  timer.registrate( 8, "   dispatch,connect_to_top");
  timer.registrate( 9, "   dispatch,symmetrize_spins");
  timer.registrate(10, "   dispatch,pack_tree");
  timer.registrate(11, "   dispatch,assign_cluster_id");
  timer.registrate(12, "   dispatch,accumulate");
  timer.registrate(13, "   dispatch,collect");
  timer.registrate(14, "   dispatch,determine_flip");
  timer.registrate(15, "   dispatch,flip_operator&spins");
  timer.registrate(16, "   dispatch,measurement");
  unifier_i.init_timer(timer);

  // initialize parallel cluster unifier
  if (enable_improved_estimator)
    unifier_i.initialize(timer, num_sites(lattice.vg()), p.value_or_default("PARTITION", ""),
      p.value_or_default("DUPLEX", true), p.value_or_default("RESERVE_CHUNK_LINKS", 0),
      p.value_or_default("RESERVE_CHUNK_ESTIMATES", 0));
  else
    unifier_n.initialize(timer, num_sites(lattice.vg()), p.value_or_default("PARTITION", ""),
      p.value_or_default("DUPLEX", true), p.value_or_default("RESERVE_CHUNK_LINKS", 0),
      p.value_or_default("RESERVE_CHUNK_ESTIMATES", 0));

  // initialize estimators
  estimator.initialize(p, lattice, model.is_signed(), enable_improved_estimator);

  timer.start(1);

#ifdef STD_OUTPUT
  if (comm.rank() == 0) {
    std::cout << std::setprecision(12);
    std::cout << "initialization done\n";
  }
#endif
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
      if (enable_improved_estimator) {
        obs << alps::RealObservable("Weight of Zero-Meron Sector");
        obs << alps::RealObservable("Sign in Zero-Meron Sector");
      }
    }
    estimator.init_observables(obs, model.is_signed(), enable_improved_estimator);
  }
}

void loop_worker::run(alps::ObservableSet& obs) {
#ifdef STD_OUTPUT
  if (comm.rank() == 0) std::cout << "MCS: " << mcs() << ' ';
#endif
  timer.start(2);
  beta = 1.0 / temperature(mcs());
  tau0 = 1.0 * comm.rank() / comm.size();
  tau1 = 1.0 * (comm.rank() + 1) / comm.size();

  //       FIELD               SIGN                IMPROVE
#ifdef LOOPER_OPENMP
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs, coll_ig, estimates_ig,
    unifier_i);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs, coll_ng, estimates_ng,
    unifier_n);
#else
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs, coll_i, estimates_i,
    unifier_i);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs, coll_n, estimates_n,
    unifier_n);
#endif

  timer.stop(2);
#ifdef STD_OUTPUT
  if (comm.rank() == 0) std::cout << std::endl << std::flush;
  timer.summarize(std::cerr);
#endif
  if (!mcs.is_thermalized()) timer.clear();
  ++mcs;
}


#ifdef LOOPER_OPENMP
template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE,
  typename UNIFIER>
void loop_worker::dispatch(alps::ObservableSet& obs, std::vector<COLLECTOR>& coll_g,
  std::vector<std::vector<ESTIMATE> >& estimates_g, UNIFIER& unifier) {
#else
template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE,
  typename UNIFIER>
void loop_worker::dispatch(alps::ObservableSet& obs, COLLECTOR& coll,
  std::vector<ESTIMATE>& estimates, UNIFIER& unifier) {
#endif
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      enable_improved_estimator != IMPROVE()) return;
  typedef COLLECTOR collector_t;
  typedef ESTIMATE estimate_t;

  timer.start(3);
  int comm_rank = comm.rank(); // make temperary since comm.rank() and comm.size() are
  int comm_size = comm.size(); // not thread-safe
#ifdef LOOPER_OPENMP
  int num_threads = omp_get_max_threads();
#endif
  int nrs = num_sites(lattice.rg());
  int nvs = num_sites(lattice.vg());

  //
  // diagonal update and cluster construction
  //

  // initialize spin & operator information
  timer.start(4);
  #ifdef LOOPER_OPENMP
  looper::padded_vector<std::vector<double>, 4> current_times(num_threads); // current time of each thread
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int tid = omp_get_thread_num();
    std::vector<local_operator_t>& operators = operators_g[tid];
    std::vector<local_operator_t>& operators_p = operators_pg[tid];
    current_times[tid] = tau0;
    #endif
    std::swap(operators, operators_p); operators.resize(0);
    // insert a diagonal operator at the end of operators_p
    operators_p.push_back(local_operator_t(0, local_operator_t::location_t(), tau1));
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) spins_c[s] = spins[s];
  } // end omp parallel
  timer.stop(4);

  // fill times
  timer.start(5);
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int tid = omp_get_thread_num();
    std::vector<double>& times = times_g[tid];
    alps::rng_helper::generator_type generator = generator_01(tid);
    expdist_t expdist(beta * model.graph_weight(tid));
    #else
    alps::rng_helper::generator_type generator = generator_01();
    expdist_t expdist(beta * model.graph_weight());
    #endif
    times.resize(0);
    double t = tau0;
    while (t < tau1) {
      t += expdist(generator);
      times.push_back(t);
    } // a sentinel (t >= tau1) will be appended
  } // end omp parallel
  timer.stop(5);

  // initialize cluster information (setup cluster fragments)
  timer.start(6);
  int boundary_offset = nvs;
  int local_offset = boundary_offset + nvs - vs_local_begin;
#ifdef LOOPER_OPENMP
  fragment_offset_g[0] = boundary_offset + nvs + nvs_local;
  for (int p = 1; p < num_threads + 1; ++p) {
    int n = operators_pg[p-1].size() + times_g[p-1].size();
    fragment_offset_g[p] = fragment_offset_g[p-1] + n;
  }
  looper::expand(fragments, fragment_offset_g[num_threads]);
#else
  int fragment_offset = boundary_offset + nvs + nvs_local;
  looper::expand(fragments, fragment_offset + operators_p.size() + times.size());
#endif
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int& fragment_offset = fragment_offset_g[0];
    #endif
    cluster_fragment_t fragment_init;
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static) nowait
    #endif
    for (int s = 0; s < fragment_offset; ++s) fragments[s] = fragment_init;
    // not necessary??
    // #ifdef LOOPER_OPENMP
    // #pragma omp for schedule(static) nowait
    // #endif
    // for (int s = boundary_offset; s < boundary_offset + nvs; ++s) fragments[s] = fragment_init;
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static) nowait
    #endif
    for (int s = 0; s < vs_local_begin; ++s) current[s] = s;
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static) nowait
    #endif
    for (int s = vs_local_begin; s < vs_local_end; ++s) current[s] = local_offset + s;
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = vs_local_end; s < nvs; ++s) current[s] = s;
  } // end omp parallel
  timer.stop(6);

  timer.start(7);
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int tid = omp_get_thread_num();
    std::vector<local_operator_t>& operators = operators_g[tid];
    std::vector<local_operator_t>& operators_p = operators_pg[tid];
    collector_t coll = coll_g[tid]; // use copy instead of reference to avoid false sharing
    std::vector<double>& times = times_g[tid];
    int& num_fragments = num_fragments_g[tid];
    int& fragment_offset = fragment_offset_g[tid];
    alps::rng_helper::generator_type generator = generator_01(tid);
    #else
    alps::rng_helper::generator_type generator = generator_01();
    #endif

    // intialize measurements
    coll.reset(estimator);
    looper::normal_accumulator<estimator_t, IMPROVE> accum_n(coll, lattice, estimator);
    if (comm_rank == 0) {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) accum_n.start_bottom(time_t(tau0), s, spins_c[s]);
    } else {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) accum_n.start(time_t(tau0), s, spins_c[s]);
    }

    int fid = fragment_offset;
    std::vector<double>::iterator tmi = times.begin();
    for (operator_iterator opi = operators_p.begin(); opi != operators_p.end();) {
      // diagonal update & labeling
      if (*tmi < opi->time()) {
        #ifdef LOOPER_OPENMP
        current_times[tid] = *tmi;
        loop_graph_t g = model.choose_graph(generator, tid);
        #else
        loop_graph_t g = model.choose_graph(generator);
        #endif
        if (is_bond(g)) {
          #ifdef LOOPER_OPENMP
          // wait for other threads
          int nid = sharing(pos(g));
          if (nid != tid) {
            do {
              #pragma omp flush (current_times)
            } while (current_times[nid] < *tmi);
          }
          #endif
          if (is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                            spins_c[target(pos(g), lattice.vg())])) {
            operators.push_back(local_operator_t(g, *tmi));
            ++tmi;
          } else {
            ++tmi;
            continue;
          }
        } else {
          operators.push_back(local_operator_t(g, *tmi));
          ++tmi;
        }
      } else {
        #ifdef LOOPER_OPENMP
        current_times[tid] = opi->time();
        #endif
        if (opi->is_diagonal()) {
          ++opi;
          continue;
        } else {
          operators.push_back(*opi);
          #ifdef LOOPER_OPENMP
          if (opi->is_bond()) {
            // wait for other threads
            int nid = sharing(opi->pos());
            if (nid != tid) {
              do {
                #pragma omp flush (current_times)
              } while (current_times[nid] < opi->time());
            }
          }
          #endif
          ++opi;
        }
      }

      operator_iterator oi = operators.end() - 1;
      if (oi->is_bond()) {
        int b = oi->pos();
        int s0 = source(b, lattice.vg());
        int s1 = target(b, lattice.vg());
        if (oi->is_offdiagonal()) {
          oi->assign_graph(model.choose_offdiagonal(generator, oi->loc(),
            spins_c[s0], spins_c[s1]));
          accum_n.end_b(oi->time(), b, s0, s1, spins_c[s0], spins_c[s1]);
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
          accum_n.begin_b(oi->time(), b, s0, s1, spins_c[s0], spins_c[s1]);
        }
        boost::tie(fid, current[s0], current[s1], oi->loop0, oi->loop1) =
          reconnect(fragments, fid, oi->graph(), current[s0], current[s1]);
      } else {
        int s = oi->pos();
        if (oi->is_offdiagonal()) {
          accum_n.end_s(oi->time(), s, spins_c[s]);
          spins_c[s] ^= 1;
          accum_n.begin_s(oi->time(), s, spins_c[s]);
        }
        boost::tie(fid, current[s], oi->loop0, oi->loop1) =
          reconnect(fragments, fid, oi->graph(), current[s]);
      }
    }
    if (comm_rank == (comm_size - 1)) {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) accum_n.stop_top(time_t(tau1), s, spins_c[s]);
    } else {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) accum_n.stop(time_t(tau1), s, spins_c[s]);
    }
    num_fragments = fid - fragment_offset;
    #ifdef LOOPER_OPENMP
    coll_g[tid] = coll; // write back to global array
    #endif
  }
  timer.stop(7);

  // connect to top
  timer.start(8);
  #ifdef LOOPER_OPENMP
  #pragma omp parallel for schedule(static)
  #endif
  for (int s = 0; s < nvs; ++s) unify(fragments, current[s], boundary_offset + s);
  timer.stop(8);

  // symmetrize spins
  timer.start(9);
  if (max_virtual_sites(lattice) == 1) {
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = vs_local_begin; s < vs_local_end; ++s) unify(fragments, s, local_offset + s);
  } else {
    #ifdef LOOPER_OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef LOOPER_OPENMP
      int tid = omp_get_thread_num();
      std::vector<int>& perm = perm_g[tid];
      alps::rng_helper::generator_type generator = generator_01(tid);
      #else
      alps::rng_helper::generator_type generator = generator_01();
      #endif
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int rs = rs_local_begin; rs < rs_local_end; ++rs) {
        looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
        boost::tie(vsi, vsi_end) = sites(lattice, rs);
        int offset1 = *vsi;
        int offset2 = *vsi - vs_local_begin;
        int s2 = *vsi_end - *vsi;
        for (int i = 0; i < s2; ++i) perm[i] = i;
        looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
          spins_t.begin() + offset2, spins.begin() + offset1, generator);
        for (int i = 0; i < s2; ++i)
          unify(fragments, offset1 + i, local_offset + vs_local_begin + offset2 + perm[i]);
      }
    }
  }
  timer.stop(9);

  timer.start(10);
  if (comm_size == 1) {
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) unify(fragments, s, boundary_offset + s);
  } else {
    pack_tree(fragments, 2 * nvs);
  }
  timer.stop(10);

  //
  // cluster flip
  //

  // assign cluster id
  timer.start(11);
  int nc, noc;
#ifdef LOOPER_OPENMP
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int ncl;
    noc_g[tid] = count_root_p(fragments, 0, nvs) +
      count_root_p(fragments, boundary_offset, nvs);
    ncc_g[tid] = count_root_p(fragments, local_offset + vs_local_begin, nvs_local) +
      count_root(fragments, fragment_offset_g[tid], num_fragments_g[tid]);
    #pragma omp barrier
    ncl = looper::subaccumulate(noc_g, tid);
    ncl = set_id_p(fragments, 0, nvs, ncl);
    ncl = set_id_p(fragments, boundary_offset, nvs, ncl);
    if (tid + 1 == num_threads) noc = ncl;
    #pragma omp barrier
    ncl = noc + looper::subaccumulate(ncc_g, tid);
    ncl = set_id_p(fragments, local_offset + vs_local_begin, nvs_local, ncl);
    ncl = set_id(fragments, fragment_offset_g[tid], num_fragments_g[tid], ncl);
    if (tid + 1 == num_threads) nc = ncl;
    #pragma omp barrier
    copy_id_p(fragments, 0, nvs);
    copy_id_p(fragments, boundary_offset, nvs);
    copy_id_p(fragments, local_offset + vs_local_begin, nvs_local);
    copy_id(fragments, fragment_offset_g[tid], num_fragments_g[tid]);
  }
  if (comm_size == 1) noc = 0;
#else
  nc = set_id(fragments, 0, nvs, 0);
  nc = set_id(fragments, boundary_offset, nvs, nc);
  noc = (comm_size == 1) ? 0 : nc;
  nc = set_id(fragments, local_offset + vs_local_begin, nvs_local, nc);
  nc = set_id(fragments, fragment_offset, num_fragments, nc);
  copy_id(fragments, 0, nvs);
  copy_id(fragments, boundary_offset, nvs);
  copy_id(fragments, local_offset + vs_local_begin, nvs_local);
  copy_id(fragments, fragment_offset, num_fragments);
#endif
  timer.stop(11);

  // accumulate physical property of clusters
  timer.start(12);
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int tid = omp_get_thread_num();
    current_times[tid] = tau0;
    std::vector<local_operator_t>& operators = operators_g[tid];
    std::vector<estimate_t>& estimates = estimates_g[tid];
    #endif
  looper::expand(estimates, nc);
    for (int c = 0; c < nc; ++c) {
      estimates[c].reset(estimator);
    }
    if (IMPROVE() || FIELD()) {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) spins_c[s] = spins[s];
      looper::improved_accumulator<estimator_t, cluster_fragment_t, ESTIMATE, IMPROVE>
        accum_i(estimates, lattice, estimator, fragments);
      if (comm_rank == 0) {
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static) nowait
        #endif
        for (int s = 0; s < vs_local_begin; ++s)
          accum_i.start_bottom(s, time_t(tau0), s, spins_c[s]);
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static) nowait
        #endif
        for (int s = vs_local_begin; s < vs_local_end; ++s)
          accum_i.start_bottom(local_offset+s, time_t(tau0), s, spins_c[s]);
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int s = vs_local_end; s < nvs; ++s)
          accum_i.start_bottom(s, time_t(tau0), s, spins_c[s]);
      } else {
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static) nowait
        #endif
        for (int s = 0; s < vs_local_begin; ++s)
          accum_i.start(s, time_t(tau0), s, spins_c[s]);
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static) nowait
        #endif
        for (int s = vs_local_begin; s < vs_local_end; ++s)
          accum_i.start(local_offset+s, time_t(tau0), s, spins_c[s]);
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int s = vs_local_end; s < nvs; ++s)
          accum_i.start(s, time_t(tau0), s, spins_c[s]);
      }
      for (operator_iterator opi = operators.begin(); opi != operators.end(); ++opi) {
        time_t t = opi->time();
        #ifdef LOOPER_OPENMP
        current_times[tid] = t;
        #endif
        if (opi->is_bond()) {
          if (!opi->is_frozen_bond_graph()) {
            int b = opi->pos();
            int s0 = source(b, lattice.vg());
            int s1 = target(b, lattice.vg());
            #ifdef LOOPER_OPENMP
            // wait for other threads
            int nid = sharing(b);
            if (nid != tid) {
              do {
                #pragma omp flush (current_times)
              } while (current_times[nid] < t);
            }
            #endif
            accum_i.end_b(opi->loop_l0(), opi->loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
            if (opi->is_offdiagonal()) {
              spins_c[s0] ^= 1;
              spins_c[s1] ^= 1;
            }
            accum_i.begin_b(opi->loop_u0(), opi->loop_u1(), t, b, s0, s1, spins_c[s0],
              spins_c[s1]);
          }
        } else {
          if (!opi->is_frozen_site_graph()) {
            //// not thread safe !!!!
            int s = opi->pos();
            accum_i.end_s(opi->loop_l(), t, s, spins_c[s]);
            if (opi->is_offdiagonal()) {
              spins_c[s] ^= 1;
            }
            accum_i.begin_s(opi->loop_u(), t, s, spins_c[s]);
          }
        }
      }
      #ifdef LOOPER_OPENMP
      current_times[tid] = tau1;
      #endif
      if (comm_rank == (comm_size - 1)) {
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int s = 0; s < nvs; ++s)
          accum_i.stop_top(boundary_offset + s, time_t(tau1), s, spins_c[s]);
      } else {
        #ifdef LOOPER_OPENMP
        #pragma omp for schedule(static)
        #endif
        for (int s = 0; s < nvs; ++s)
          accum_i.stop(boundary_offset + s, time_t(tau1), s, spins_c[s]);
      }
    }
  }
  timer.stop(12);

  // accumulate cluster properties
  timer.start(13);
#ifdef LOOPER_OPENMP
  if (IMPROVE()) {
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      collector_t coll = coll_g[tid]; // use copy instead of reference to avoid false sharing
      coll.set_num_operators(operators_g[tid].size());
      #pragma omp for schedule(static)
      for (int c = 0; c < nc; ++c) {
        for (int p = 1; p < num_threads; ++p) estimates_g[0][c] += estimates_g[p][c];
      }
      #pragma omp for schedule(static)
      for (int c = noc; c < nc; ++c) coll += estimates_g[0][c];
      coll_g[tid] = coll; // write back to global array
    } // end omp parallel
  } else {
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      coll_g[tid].set_num_operators(operators_g[tid].size());
    }
  }
  std::vector<estimate_t>& estimates = estimates_g[0];
  collector_t& coll = coll_g[0];
  for (int p = 1; p < num_threads; ++p) coll += coll_g[p];
#else
  coll.set_num_operators(operators.size());
  if (IMPROVE()) for (int c = noc; c < nc; ++c) coll += estimates[c];
#endif
  coll.set_num_open_clusters(noc);
  coll.set_num_clusters(nc - noc);
  timer.stop(13);

  // determine whether clusters are flipped or not
  timer.start(14);
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int tid = omp_get_thread_num();
    alps::rng_helper::generator_type generator = generator_01(tid);
    #else
    alps::rng_helper::generator_type generator = generator_01();
    #endif
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int c = 0; c < nc; ++c) estimates[c].to_flip = (generator() < 0.5);
  }
  if (comm_size > 1) {
    #ifdef LOOPER_OPENMP
    unifier.unify(coll, fragments, 0, boundary_offset, estimates_g, timer);
    #else
    unifier.unify(coll, fragments, 0, boundary_offset, estimates, timer);
    #endif
  }
  timer.stop(14);

  // flip operators & spins
  timer.start(15);
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef LOOPER_OPENMP
    int tid = omp_get_thread_num();
    std::vector<local_operator_t>& operators = operators_g[tid];
    #endif
    for (operator_iterator opi = operators.begin(); opi != operators.end(); ++opi) {
      if ((estimates[fragments[opi->loop_0()].id()].to_flip ^
           estimates[fragments[opi->loop_1()].id()].to_flip) & 1) opi->flip();
    }
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = 0; s < vs_local_begin; ++s)
      if (estimates[fragments[s].id()].to_flip & 1) spins[s] ^= 1;
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = vs_local_begin; s < vs_local_end; ++s) {
      if (estimates[fragments[s].id()].to_flip & 1) spins_t[s - vs_local_begin] ^= 1;
      if (estimates[fragments[local_offset + s].id()].to_flip & 1) spins[s] ^= 1;
    }
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = vs_local_end; s < nvs; ++s)
      if (estimates[fragments[s].id()].to_flip & 1) spins[s] ^= 1;
  }
  timer.stop(15);

  //
  // measurement
  //

  timer.start(16);
  if (comm_rank == 0) {
    obs["Temperature"] << 1/beta;
    obs["Inverse Temperature"] << beta;
    obs["Volume"] << (double)lattice.volume();
    obs["Number of Sites"] << (double)nrs;
    obs["Number of Clusters"] << coll.num_clusters();
#ifdef STD_OUTPUT
    std::cout << 1/beta << ' '
              << beta << ' '
              << (double)lattice.volume() << ' '
              << (double)num_sites(lattice.rg()) << ' '
              << coll.num_clusters() << ' ';
#endif
    double nop = coll.num_operators();
    double ene = model.energy_offset() - nop / beta;
    coll.set_energy(ene);

    coll.commit(obs, estimator, lattice, beta, 1, nop, spins);
  }
  timer.stop(16);
  timer.stop(3);
}

boost::tuple<int, int> loop_worker::real_site_range(int nrs, int nprocs, int rank) const {
  int nrs_local = nrs / nprocs;
  int nextra = nrs % nprocs;
  int rs_begin = rank * nrs_local + (rank < nextra ? rank : nextra);
  int rs_local_end = rs_begin + nrs_local + (rank < nextra ? 1 : 0);
  return boost::make_tuple(rs_begin, rs_local_end);
}

typedef looper::evaluator<loop_config::measurement_set> loop_evaluator;

//
// dynamic registration to the factories
//

PARAPACK_REGISTER_PARALLEL_ALGORITHM(loop_worker, "loop; path integral");

} // end namespace
