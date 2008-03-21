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

#include <parapack/worker.h>

inline int log2(int n) {
  int r = -1;
  for (; n > 0; n = (n >> 1)) ++r;
  return r;
}

struct walker_direc {
  enum walker_direc_t { up, down, unlabeled };
};
typedef walker_direc::walker_direc_t walker_direc_t;

class loop_worker :
  public alps::parapack::mc_worker, public alps::parapack::need_multiple_observableset_tag {
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

  loop_worker(alps::Parameters const& p, std::vector<alps::ObservableSet>& obs);
  virtual ~loop_worker() {}

  static void print_copyright(std::ostream&) {}

  void run(std::vector<alps::ObservableSet>& obs);

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void save(alps::ODump& dp) const {
    uint32_t direc_tmp = (direc == walker_direc::up ? 0 : 1);
    dp << mcs << spins << operators << logg << direc_tmp << factor;
  }
  void load(alps::IDump& dp) {
    uint32_t direc_tmp;
    dp >> mcs >> spins >> operators >> logg >> direc_tmp >> factor;
    direc = (direc_tmp == 0 ? walker_direc::up : walker_direc::down);
  }

protected:
  void build(std::vector<alps::ObservableSet>& obs);

  template<typename FIELD, typename SIGN, typename IMPROVE>
  void flip(std::vector<alps::ObservableSet>& obs);

  void measure(std::vector<alps::ObservableSet>& obs);

private:
  // helpers
  lattice_t lattice;
  model_t model;

  // parameters
  bool use_improved_estimator;
  bool conventional;
  int max_order;
  int interval;
  double final, flatness;
  bool correct;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;
  std::valarray<double> logg;
  walker_direc_t direc;
  double factor;

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

loop_worker::loop_worker(alps::Parameters const& p, std::vector<alps::ObservableSet>& obs) :
  alps::parapack::mc_worker(p),
  lattice(p), model(p, lattice, looper::is_path_integral<mc_type>::type()), mcs(p) {

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

  conventional = p.value_or_default("CONVENTIONAL_QUANTUM_WANG_LANDAU", true);
  max_order = static_cast<int>(evaluate("MAX_EXPANSION_ORDER", p));
  factor = p.defined("INITIAL_UPDATE_FACTOR") ? std::log(evaluate("INITIAL_UPDATE_FACTOR", p)) : 1;
  if (conventional) {
    interval = p.defined("CHECK_INTERVAL") ? static_cast<int>(evaluate("CHECK_INTERVAL", p)) : 1024;
    final = std::log(static_cast<double>(p.value_or_default("FINAL_UPDATE_FACTOR", 1.00000001)));
    mcs.set_thermalization(interval);
    flatness = p.value_or_default("FLATNESS_THRESHOLD", 0.95);
  }
  correct = p.value_or_default("FACTOR_CORRECTION", false);

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);
  logg.resize(max_order);
  for (int i = 0; i < max_order; ++i) logg[i] = -2 * i;
  direc = walker_direc::down;

  int num_part = log2(max_order-1) + 2;
  obs.resize(num_part);
  obs[0] << alps::HistogramObservable<int>("Global Histogram", 0, max_order);
  obs[0] << alps::HistogramObservable<int>("Histogram of Upward-moving Walker", 0, max_order);
  obs[0] << alps::HistogramObservable<int>("Local Histogram", 0, max_order);
  obs[0] << alps::HistogramObservable<int>("Partition Histogram", 0, num_part);
  obs[0] << make_observable(alps::RealObservable("Inverse Round-trip Time"));
  obs[0] << alps::SimpleRealVectorObservable("Log(g)");

  BOOST_FOREACH(alps::ObservableSet& o, obs) {
    o << make_observable(alps::SimpleRealObservable("Volume"));
    o << make_observable(alps::SimpleRealObservable("Number of Sites"));
    o << make_observable(alps::SimpleRealObservable("Number of Clusters"));
    if (model.is_signed()) {
      o << alps::RealObservable("Sign");
      if (use_improved_estimator) {
        o << alps::RealObservable("Weight of Zero-Meron Sector");
        o << alps::RealObservable("Sign in Zero-Meron Sector");
      }
    }
    estimator.initialize(o, p, lattice, model.is_signed(), use_improved_estimator);
  }
}

void loop_worker::run(std::vector<alps::ObservableSet>& obs) {
  ++mcs;

  build(obs);

  //   FIELD               SIGN                IMPROVE
  flip<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::true_ >(obs);
  flip<boost::mpl::true_,  boost::mpl::true_,  boost::mpl::false_>(obs);
  flip<boost::mpl::true_,  boost::mpl::false_, boost::mpl::true_ >(obs);
  flip<boost::mpl::true_,  boost::mpl::false_, boost::mpl::false_>(obs);
  flip<boost::mpl::false_, boost::mpl::true_,  boost::mpl::true_ >(obs);
  flip<boost::mpl::false_, boost::mpl::true_,  boost::mpl::false_>(obs);
  flip<boost::mpl::false_, boost::mpl::false_, boost::mpl::true_ >(obs);
  flip<boost::mpl::false_, boost::mpl::false_, boost::mpl::false_>(obs);

  measure(obs);

  if (conventional && mcs() == mcs.thermalization()) {
    alps::HistogramObservable<int> *hist =
      dynamic_cast<alps::HistogramObservable<int>*>(&obs[0]["Local Histogram"]);
    double mean = 0;
    unsigned int hist_min = mcs();
    for (int i = 0; i < max_order; ++i) {
      mean += (*hist)[i];
      hist_min = std::min(hist_min, (*hist)[i]);
    }
    mean /= max_order;
    if (hist_min > flatness * mean) {
      std::clog << mcs() << ": flatness check succeeded (mean = " << mean << ", Hmin = "
                << hist_min << ", ratio = " << hist_min / mean << ")\n";
      factor *= 0.5;
      if (factor > final) {
        std::clog << mcs() << ": update factor is reduced to exp(" << factor
                  << ") (target = exp(" << final << "))\n";
        mcs.set_thermalization(mcs.thermalization() + interval);
        hist->reset(false);
      } else {
        std::clog << mcs() << ": thermalization done\n";
        factor = 0;
      }
    } else {
      std::clog << mcs() << ": flatness check failed (mean = " << mean << ", Hmin = "
                << hist_min << ", ratio = " << hist_min / mean << ")\n";
      mcs.set_thermalization(mcs.thermalization() + interval);
    }
  }
}


//
// diagonal update and cluster construction
//

void loop_worker::build(std::vector<alps::ObservableSet>& obs) {

  alps::HistogramObservable<int> *hist =
    dynamic_cast<alps::HistogramObservable<int>*>(&obs[0]["Local Histogram"]);

  // initialize spin & operator information
  int nop = operators.size();
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double bw = model.graph_weight();
  double factor2 = correct ? 0.5 * factor : 0;
  bool try_gap = true;
  for (operator_iterator opi = operators_p.begin(); try_gap || opi != operators_p.end();) {

    bool has_operator;

    // diagonal update & labeling
    if (try_gap) {
      has_operator = false;
      if (nop < max_order-1 &&
          (nop + 1) * uniform_01() < bw * std::exp(logg[nop] - (logg[nop+1] + factor2))) {
        loop_graph_t g = model.choose_graph(uniform_01);
        if ((is_bond(g) && is_compatible(g, spins_c[source(pos(g), lattice.vg())],
                                            spins_c[target(pos(g), lattice.vg())])) ||
            (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
          operators.push_back(local_operator_t(g));
          ++nop;
          has_operator = true;
        }
      }
      if (!has_operator) try_gap = false;
    } else {
      has_operator = true;
      if (opi->is_diagonal()) {
        if (bw * uniform_01() < nop * std::exp(logg[nop] - (logg[nop-1] + factor2))) {
          --nop;
          ++opi;
          has_operator = false;
        } else {
          if (opi->is_site()) {
            opi->assign_graph(model.choose_diagonal(uniform_01, opi->loc(), spins_c[opi->pos()]));
          } else {
            opi->assign_graph(model.choose_diagonal(uniform_01, opi->loc(),
              spins_c[source(opi->pos(), lattice.vg())],
              spins_c[target(opi->pos(), lattice.vg())]));
          }
        }
      } else {
        if (opi->is_bond())
          opi->assign_graph(model.choose_offdiagonal(uniform_01, opi->loc(),
            spins_c[source(opi->pos(), lattice.vg())],
            spins_c[target(opi->pos(), lattice.vg())]));
      }
      if (has_operator) {
        operators.push_back(*opi);
        ++opi;
        try_gap = true;
      }
    }
    //// std::cerr << "visit " << nop << std::endl;
    if (hist) *hist << nop;
    logg[nop] += factor;
    if (direc == walker_direc::up)
      obs[0]["Histogram of Upward-moving Walker"] << nop;
    if (direc == walker_direc::up && nop == 0) {
      obs[0]["Inverse Round-trip Time"] << 1.;
      direc = walker_direc::down;
    } else {
      obs[0]["Inverse Round-trip Time"] << 0.;
    }
    if (direc == walker_direc::down && nop == max_order - 1)
      direc = walker_direc::up;

    if (!has_operator) continue;

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

  // normalize extended weight
  double logg0 = logg[0];
  for (int i = 0; i < max_order; ++i) logg[i] -= logg0;

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
        spins.begin() + offset, spins_c.begin() + offset, uniform_01);
      for (int i = 0; i < s2; ++i) unify(fragments, offset+i, current[offset+perm[i]]);
    }
  }
}


//
// cluster flip
//

template<typename FIELD, typename SIGN, typename IMPROVE>
void loop_worker::flip(std::vector<alps::ObservableSet>& obs) {
  if (model.has_field() != FIELD() ||
      model.is_signed() != SIGN() ||
      use_improved_estimator != IMPROVE()) return;

  int nvs = num_sites(lattice.vg());
  int nop = operators.size();
  int part = log2(nop) + 1;

  // assign cluster id
  int nc = 0;
  BOOST_FOREACH(cluster_fragment_t& f, fragments) if (f.is_root()) f.set_id(nc++);
  BOOST_FOREACH(cluster_fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
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
    ci.to_flip = ((2*uniform_01()-1) < (FIELD() ? std::tanh(ci.weight) : 0));
    if (SIGN() && IMPROVE() && (ci.sign & 1 == 1)) improved_sign = 0;
  }

  // improved measurement
  if (IMPROVE()) {
    typename looper::collector<estimator_t>::type coll = get_collector(estimator);
    coll = std::accumulate(estimates.begin(), estimates.end(), coll);
    estimator.improved_measurement(obs[part], lattice, 1, improved_sign, spins, operators,
      spins_c, fragments, coll);
    if (SIGN()) {
      obs[part]["Sign"] << improved_sign;
      if (alps::is_zero(improved_sign)) {
        obs[part]["Weight of Zero-Meron Sector"] << 0.;
      } else {
        obs[part]["Weight of Zero-Meron Sector"] << 1.;
        obs[part]["Sign in Zero-Meron Sector"] << improved_sign;
      }
    }
  }
  obs[part]["Number of Clusters"] << (double)clusters.size();

  // flip operators & spins
  BOOST_FOREACH(local_operator_t& op, operators)
    if (clusters[fragments[op.loop_0()].id()].to_flip ^
        clusters[fragments[op.loop_1()].id()].to_flip)
      op.flip();
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id()].to_flip) spins[s] ^= 1;
}


//
// measurement
//

void loop_worker::measure(std::vector<alps::ObservableSet>& obs) {

  int nop = operators.size();
  int part = log2(nop) + 1;

  obs[part]["Volume"] << (double)lattice.volume();
  obs[part]["Number of Sites"] << (double)num_sites(lattice.rg());

  obs[0]["Global Histogram"] << nop;
  obs[0]["Partition Histogram"] << part;
  if (progress() >= 1) obs[0]["Log(g)"] << logg;

  // sign
  if (model.is_signed() && !use_improved_estimator) obs[part]["Sign"] << sign;

  // other quantities
  estimator.normal_measurement(obs[part], lattice, 1, sign, spins, operators, spins_c);
}
