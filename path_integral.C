/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_config.h"
#include <looper/capacity.h>
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/expand.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/temperature.h>
#include <looper/timer.hpp>
#include <looper/type.h>
#include <looper/union_find.h>
#include <alps/parapack/worker.h>
#include <alps/parapack/exchange.h>
#include <alps/numeric/is_zero.hpp>

#ifdef LOOPER_OPENMP
# include <looper/padded_vector.h>
# include <looper/subaccumulate.h>
# include <omp.h>
#endif

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
  typedef double weight_parameter_type;

  typedef boost::exponential_distribution<> expdist_t;

  loop_worker(alps::Parameters const& p);
  virtual ~loop_worker() {
#ifdef LOOPER_OPENMP
    if (reserved) {
      looper::vector_capacity capacity(times_g, operators_g, operators_pg, estimates_ig, estimates_ng, fragments);
      capacity.report();
    }
#else
    if (reserved) {
      looper::vector_capacity capacity(times, operators, operators_p, estimates_i, estimates_n, fragments);
      capacity.report();
    }
#endif
    timer.stop(1);
    timer.summarize();
  }

  void init_observables(alps::Parameters const& params, alps::ObservableSet& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(alps::ObservableSet& obs);

  // for exchange Monte Carlo
  void set_beta(double beta) { temperature.set_beta(beta); }
  double weight_parameter() const {
#ifdef LOOPER_OPENMP
    int num_threads = omp_get_max_threads();
    int n = 0;
    for (int p = 0; p < num_threads; ++p) n += operators_g[p].size();
    return n;
#else
    return operators.size();
#endif
  }
  static double log_weight(double gw, double beta) { return std::log(beta) * gw; }

  void save(alps::ODump& dp) const {
#ifdef LOOPER_OPENMP
    dp << mcs << spins << operators_g;
#else
    dp << mcs << spins << operators;
#endif
  }
  void load(alps::IDump& dp) {
#ifdef LOOPER_OPENMP
    dp >> mcs >> spins >> operators_g;
#else
    dp >> mcs >> spins >> operators;
#endif
  }

protected:
#ifdef LOOPER_OPENMP
  template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
  void dispatch(alps::ObservableSet& obs, std::vector<COLLECTOR>& coll,
    std::vector<std::vector<ESTIMATE> >& estimates);
#else
  template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
  void dispatch(alps::ObservableSet& obs, COLLECTOR& coll, std::vector<ESTIMATE>& estimates);
#endif

private:
  // helpers
  lattice_t lattice;
#ifdef LOOPER_OPENMP
  looper::simple_lattice_sharing sharing;
#endif
  model_t model;

  // parameters
  looper::temperature temperature;
  double beta;
  bool enable_improved_estimator;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
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
  std::vector<int> nc_g;
  std::vector<std::vector<cluster_info_t> > clusters_g;
  std::vector<std::vector<estimator_t::improved_estimator::estimate> > estimates_ig;
  std::vector<std::vector<estimator_t::normal_estimator::estimate> > estimates_ng;
  std::vector<std::vector<int> > perm_g;
#else
  std::vector<local_operator_t> operators_p;
  std::vector<double> times;
  int fragment_offset;
  int num_fragments;
  std::vector<cluster_info_t> clusters;
  std::vector<estimator_t::improved_estimator::estimate> estimates_i;
  std::vector<estimator_t::normal_estimator::estimate> estimates_n;
  std::vector<int> perm;
#endif

  alps::parapack::timer timer;
  bool reserved;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::Parameters const& p)
  : alps::parapack::mc_worker(p), lattice(p),
#ifdef LOOPER_OPENMP
    sharing(/* will be initialized in model */),
    model(p, lattice, sharing, /* is_path_integral = */ true),
#else
    model(p, lattice, /* is_path_integral = */ true),
#endif
    temperature(p), mcs(p), timer() {

  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument("longer annealing steps than thermalization"));

  model.check_parameter(support_longitudinal_field, support_negative_sign);

  enable_improved_estimator = (!model.has_field()) && (!p.defined("DISABLE_IMPROVED_ESTIMATOR"));
  if (!enable_improved_estimator) std::cout << "WARNING: improved estimator is disabled\n";

  // configuration
#ifdef LOOPER_OPENMP
  int num_threads = omp_get_max_threads();
#endif
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
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
  fragment_offset_g.resize(num_threads + 1);
  num_fragments_g.resize(num_threads);
  nc_g.resize(num_threads);
  times_g.resize(num_threads);
  clusters_g.resize(num_threads);
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
    capacity.report();
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
  capacity.report();
#endif

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

  // initialize estimators
  estimator.initialize(p, lattice, model.is_signed(), enable_improved_estimator);

  timer.start(1);
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
  timer.start(2);
  beta = 1.0 / temperature(mcs());

  //       FIELD               SIGN                IMPROVE
#ifdef LOOPER_OPENMP
  dispatch<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::true_ >(obs, coll_ig, estimates_ig);
  dispatch<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::false_>(obs, coll_ng, estimates_ng);
  dispatch<boost::mpl::true_,  boost::mpl::false_, boost::mpl::true_ >(obs, coll_ig, estimates_ig);
  dispatch<boost::mpl::true_,  boost::mpl::false_, boost::mpl::false_>(obs, coll_ng, estimates_ng);
  dispatch<boost::mpl::false_, boost::mpl::true_,  boost::mpl::true_ >(obs, coll_ig, estimates_ig);
  dispatch<boost::mpl::false_, boost::mpl::true_,  boost::mpl::false_>(obs, coll_ng, estimates_ng);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs, coll_ig, estimates_ig);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs, coll_ng, estimates_ng);
#else
  dispatch<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::false_>(obs, coll_n, estimates_n);
  dispatch<boost::mpl::true_,  boost::mpl::false_, boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::true_,  boost::mpl::false_, boost::mpl::false_>(obs, coll_n, estimates_n);
  dispatch<boost::mpl::false_, boost::mpl::true_,  boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::false_, boost::mpl::true_,  boost::mpl::false_>(obs, coll_n, estimates_n);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs, coll_i, estimates_i);
  dispatch<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs, coll_n, estimates_n);
#endif

  timer.stop(2);
  if (!mcs.is_thermalized()) timer.clear();
  ++mcs;
}


#ifdef LOOPER_OPENMP
template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
void loop_worker::dispatch(alps::ObservableSet& obs, std::vector<COLLECTOR>& coll_g,
  std::vector<std::vector<ESTIMATE> >& estimates_g) {
#else
template<typename FIELD, typename SIGN, typename IMPROVE, typename COLLECTOR, typename ESTIMATE>
void loop_worker::dispatch(alps::ObservableSet& obs, COLLECTOR& coll,
  std::vector<ESTIMATE>& estimates) {
#endif
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      enable_improved_estimator != IMPROVE()) return;
  typedef COLLECTOR collector_t;
  typedef ESTIMATE estimate_t;

  timer.start(3);
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
    current_times[tid] = 0;
    #endif
    std::swap(operators, operators_p); operators.resize(0);
    // insert a diagonal operator at the end of operators_p
    operators_p.push_back(local_operator_t(0, local_operator_t::location_t(), 1));
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) spins_c[s] = spins[s];
  }
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
    double t = 0;
    while (t < 1) {
      t += expdist(generator);
      times.push_back(t);
    } // a sentinel (t >= 1) will be appended
  }
  timer.stop(5);

  // initialize cluster information (setup cluster fragments)
  timer.start(6);
#ifdef LOOPER_OPENMP
  fragment_offset_g[0] = nvs;
  for (int p = 1; p < num_threads + 1; ++p) {
    int n = operators_pg[p-1].size() + times_g[p-1].size();
    fragment_offset_g[p] = fragment_offset_g[p-1] + n;
  }
  looper::expand(fragments, fragment_offset_g[num_threads]);
#else
  fragment_offset = nvs;
  looper::expand(fragments, fragment_offset + operators_p.size() + times.size());
#endif
  #ifdef LOOPER_OPENMP
  #pragma omp parallel
  #endif
  {
    cluster_fragment_t fragment_init;
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) {
      fragments[s] = fragment_init;
      current[s] = s;
    }
  }
  timer.stop(6);

  timer.start(7);
  int negop = 0; // number of operators with negative weights
  #ifdef LOOPER_OPENMP
  #pragma omp parallel reduction(+:negop)
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
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) accum_n.start_bottom(time_t(0), s, spins_c[s]);

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
    #ifdef LOOPER_OPENMP
    #pragma omp for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) accum_n.stop_top(time_t(1), s, spins_c[s]);
    num_fragments = fid - fragment_offset;
    #ifdef LOOPER_OPENMP
    coll_g[tid] = coll; // write back to global array
    #endif
  }
  double sign = ((negop & 1) == 1) ? -1 : 1;
  timer.stop(7);

  timer.start(8);
  timer.stop(8);

  // symmetrize spins
  timer.start(9);
  if (max_virtual_sites(lattice) == 1) {
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int s = 0; s < nvs; ++s) unify(fragments, s, current[s]);
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
      for (int rs = 0; rs < nrs; ++rs) {
        looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
        boost::tie(vsi, vsi_end) = sites(lattice, rs);
        int offset = *vsi;
        int s2 = *vsi_end - *vsi;
        for (int i = 0; i < s2; ++i) perm[i] = i;
        looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
          spins.begin() + offset, spins_c.begin() + offset, generator);
        for (int i = 0; i < s2; ++i) unify(fragments, offset+i, current[offset+perm[i]]);
      }
    }
  }
  timer.stop(9);

  timer.start(10);
  timer.stop(10);

  //
  // cluster flip
  //

  // assign cluster id
  timer.start(11);
  int nc;
#ifdef LOOPER_OPENMP
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int ncl;
    nc_g[tid] = count_root_p(fragments, 0, nvs) +
      count_root(fragments, fragment_offset_g[tid], num_fragments_g[tid]);
    #pragma omp barrier
    ncl = looper::subaccumulate(nc_g, tid);
    ncl = set_id_p(fragments, 0, nvs, ncl);
    ncl = set_id(fragments, fragment_offset_g[tid], num_fragments_g[tid], ncl);
    if (tid + 1 == num_threads) nc = ncl;
    #pragma omp barrier
    copy_id_p(fragments, 0, nvs);
    copy_id(fragments, fragment_offset_g[tid], num_fragments_g[tid]);
  }
#else
  nc = set_id(fragments, 0, fragment_offset + num_fragments, 0);
  copy_id(fragments, 0, fragment_offset + num_fragments);
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
    current_times[tid] = 0;
    std::vector<local_operator_t>& operators = operators_g[tid];
    std::vector<cluster_info_t >& clusters = clusters_g[tid];
    std::vector<estimate_t>& estimates = estimates_g[tid];
    #endif
    looper::expand(clusters, nc);
    looper::expand(estimates, nc);
    cluster_info_t cluster_init;
    for (int c = 0; c < nc; ++c) {
      clusters[c] = cluster_init;
      estimates[c].reset(estimator);
    }
    if (IMPROVE() || FIELD()) {
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) spins_c[s] = spins[s];
      cluster_info_t::accumulator<cluster_fragment_t, FIELD, SIGN, IMPROVE>
        weight(clusters, fragments, model.field(), model.bond_sign(), model.site_sign());
      looper::improved_accumulator<estimator_t, cluster_fragment_t, ESTIMATE, IMPROVE>
        accum_i(estimates, lattice, estimator, fragments);
      #ifdef LOOPER_OPENMP
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) {
        weight.start_bottom(s, time_t(0), s, spins_c[s]);
        accum_i.start_bottom(s, time_t(0), s, spins_c[s]);
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
            weight.end_b(opi->loop_l0(), opi->loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
            accum_i.end_b(opi->loop_l0(), opi->loop_l1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
            if (opi->is_offdiagonal()) {
              spins_c[s0] ^= 1;
              spins_c[s1] ^= 1;
            }
            weight.begin_b(opi->loop_u0(), opi->loop_u1(), t, b, s0, s1, spins_c[s0], spins_c[s1]);
            accum_i.begin_b(opi->loop_u0(), opi->loop_u1(), t, b, s0, s1, spins_c[s0],
              spins_c[s1]);
          }
        } else {
          if (!opi->is_frozen_site_graph()) {
            //// not thread safe !!!!
            int s = opi->pos();
            weight.end_s(opi->loop_l(), t, s, spins_c[s]);
            accum_i.end_s(opi->loop_l(), t, s, spins_c[s]);
            if (opi->is_offdiagonal()) spins_c[s] ^= 1;
            weight.begin_s(opi->loop_u(), t, s, spins_c[s]);
            accum_i.begin_s(opi->loop_u(), t, s, spins_c[s]);
          }
        }
      }
      #ifdef LOOPER_OPENMP
      current_times[tid] = 1;
      #pragma omp for schedule(static)
      #endif
      for (int s = 0; s < nvs; ++s) {
        weight.stop_top(current[s], time_t(1), s, spins_c[s]);
        accum_i.stop_top(current[s], time_t(1), s, spins_c[s]);
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
        for (int p = 1; p < num_threads; ++p) {
          clusters_g[0][c] += clusters_g[p][c];
          estimates_g[0][c] += estimates_g[p][c];
        }
        coll += estimates_g[0][c];
      }
      coll_g[tid] = coll; // write back to global array
    } // end omp parallel
  } else {
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      coll_g[tid].set_num_operators(operators_g[tid].size());
      #pragma omp for schedule(static)
      for (int c = 0; c < nc; ++c) {
        for (int p = 1; p < num_threads; ++p) clusters_g[0][c] += clusters_g[p][c];
      }
    }
  }
  std::vector<cluster_info_t>& clusters = clusters_g[0];
  std::vector<estimate_t>& estimates = estimates_g[0];
  collector_t& coll = coll_g[0];
  for (int p = 1; p < num_threads; ++p) coll += coll_g[p];
#else
  coll.set_num_operators(operators.size());
  if (IMPROVE()) for (int c = 0; c < nc; ++c) coll += estimates[c];
#endif
  coll.set_num_clusters(nc);
  timer.stop(13);

  // determine whether clusters are flipped or not
  timer.start(14);
  negop = 0;
  #ifdef LOOPER_OPENMP
  #pragma omp parallel reduction(+:negop)
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
    for (int c = 0; c < nc; ++c) {
      estimates[c].to_flip =
        ((2*generator()-1) < (FIELD() ? std::tanh(beta * clusters[c].weight) : 0));
      if (SIGN() && IMPROVE()) negop += clusters[c].sign;
    }
  }
  double improved_sign = ((negop & 1) == 1) ? 0 : sign;
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
    for (int s = 0; s < nvs; ++s)
      if (estimates[fragments[s].id()].to_flip & 1) spins[s] ^= 1;
  }
  timer.stop(15);

  //
  // measurement
  //

  timer.start(16);
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
  if (FIELD()) {
    #ifdef LOOPER_OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (int c = 0; c < nc; ++c) {
      ene += ((estimates[c].to_flip & 1) ? -clusters[c].weight : clusters[c].weight);
    }
  }
  coll.set_energy(ene);
  coll.commit(obs, estimator, lattice, beta, improved_sign, nop, spins);
  timer.stop(16);
  timer.stop(3);
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
