/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_factory.h"
#include "loop_worker.h"
#include <looper/cluster.h>
#include <looper/histogram.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/type.h>
#include <alps/fixed_capacity_vector.h>
#include <alps/plot.h>
#include <boost/regex.hpp>

namespace {

class loop_worker : public loop_worker_base
{
public:
  typedef looper::sse                   qmc_type;
  typedef loop_worker_base              super_type;

  typedef looper::local_operator<qmc_type, loop_graph_t, time_t>
                                        local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef operator_string_t::iterator   operator_iterator;

  typedef looper::union_find::node      cluster_fragment_t;
  typedef looper::cluster_info          cluster_info_t;

  typedef loop_config::estimator_t      estimator_t;
  typedef looper::measurement::estimate<estimator_t>::type estimate_t;

  loop_worker(alps::ProcessList const& w, alps::Parameters const& p, int n);
  void dostep();

  bool is_thermalized() const { return true; }
  double work_done() const { return mcs.progress(); }

  void save(alps::ODump& dp) const
  {
    super_type::save(dp);
    dp << mcs << spins << operators << logf << histogram
       << measurement_histograms;
  }
  void load(alps::IDump& dp)
  {
    super_type::load(dp);
    dp >> mcs >> spins >> operators >> logf >> histogram
       >> measurement_histograms;
  }

protected:
  void build();
  template<typename BIPARTITE, typename FIELD, typename SIGN, typename IMPROVE>
  void flip();
  template<typename BIPARTITE, typename IMPROVE>
  void measure();

private:
  // parmeters
  looper::integer_range<int> exp_range;
  bool store_all_histograms;
  int min_visit;
  double flatness;

  looper::wl_steps mcs;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;
  double logf;
  looper::wl_histogram histogram;
  looper::histogram_set<double> measurement_histograms;

  // working vectors (no checkpointing)
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<cluster_info_t> clusters;
  std::vector<estimate_t> estimates;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::ProcessList const& w,
                         alps::Parameters const& p, int n)
  : super_type(w, p, n, looper::is_path_integral<qmc_type>::type()),
    exp_range(p.value_or_default("EXPANSION_RANGE", "[0:500]")),
    mcs(p, exp_range),
    histogram(exp_range)
{
  if (has_field())
    boost::throw_exception(std::logic_error("longitudinal field is not "
      "supported in SSE representation"));

  //
  // initialize configuration
  //

  int nvs = num_sites(vlattice());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  operators.resize(0);
  spins_c.resize(nvs);
  current.resize(nvs);

  //
  // Wang Landau parameters
  //

  if (exp_range.min() < 0)
    boost::throw_exception(std::invalid_argument("minimum of expansion order "
                                                 "must not be negative"));
  double f =
    p.value_or_default("INITIAL_MODIFICATION_FACTOR",
                       mcs.use_zhou_bhatt() ? std::exp(1) :
                       std::exp(exp_range.max() * std::log(1.*nvs) /
                                mcs.mcs_block()));
  logf = std::log(f);
  if (logf <= 0)
    boost::throw_exception(std::invalid_argument("initial modification factor "
                                                 "must be larger than 1"));
  if (mcs.use_zhou_bhatt()) {
    min_visit = static_cast<int>(1 / logf);
    flatness = parms.value_or_default("FLATNESS_THRESHOLD", -1.);
  } else {
    min_visit = 0;
    flatness = parms.value_or_default("FLATNESS_THRESHOLD", 0.2);
  }

  //
  // initialize measurements
  //

  store_all_histograms =  parms.defined("STORE_ALL_HISTOGRAMS");
  measurement_histograms.initialize(exp_range);
  measurements
    << alps::SimpleRealObservable("Energy Offset")
    << alps::SimpleRealVectorObservable("Partition Function Coefficient")
    << alps::SimpleRealVectorObservable("Histogram");
  if (store_all_histograms) {
    for (int p = 0; p < mcs.num_iterations(); ++p) {
      std::string suffix =
        "(iteration #" + boost::lexical_cast<std::string>(p) + ")";
      measurements
        << alps::SimpleRealVectorObservable("Partition Function Coefficient " +
                                            suffix)
        << alps::SimpleRealVectorObservable("Histogram " + suffix);
    }
  }
  if (is_signed()) measurement_histograms.add_histogram("Sign");
  estimator_t::initialize(measurement_histograms, is_bipartite(), is_signed(),
                          use_improved_estimator());
  measurements.reset(true);
}

void loop_worker::dostep()
{
  namespace mpl = boost::mpl;

  if (!mcs.can_work()) return;
  ++mcs;
  super_type::dostep();
  measurements["Energy Offset"] << energy_offset();

  build();

  //   BIPARTITE    FIELD        SIGN         IMPROVE
  flip<mpl::true_,  mpl::false_, mpl::true_,  mpl::true_ >();
  flip<mpl::true_,  mpl::false_, mpl::true_,  mpl::false_>();
  flip<mpl::true_,  mpl::false_, mpl::false_, mpl::true_ >();
  flip<mpl::true_,  mpl::false_, mpl::false_, mpl::false_>();
  flip<mpl::false_, mpl::false_, mpl::true_,  mpl::true_ >();
  flip<mpl::false_, mpl::false_, mpl::true_,  mpl::false_>();
  flip<mpl::false_, mpl::false_, mpl::false_, mpl::true_ >();
  flip<mpl::false_, mpl::false_, mpl::false_, mpl::false_>();

  //      BIPARTITE    IMPROVE
  measure<mpl::true_,  mpl::true_ >();
  measure<mpl::true_,  mpl::false_>();
  measure<mpl::false_, mpl::true_ >();
  measure<mpl::false_, mpl::false_>();

  if (!mcs.doing_multicanonical() && mcs() == mcs.mcs_block()) {
    if (histogram.check_flatness(flatness) &&
        histogram.check_visit(min_visit)) {
      std::cerr << "stage " << mcs.stage() << ": histogram becomes flat\n";
      histogram.subtract();
      std::string suffix =
        "(iteration #" + boost::lexical_cast<std::string>(mcs.stage()) + ")";
      histogram.store(measurements, "Partition Function Coefficient " + suffix,
                      "Histogram " + suffix, mcs.doing_multicanonical());
      logf = 0.5 * logf;
      if (mcs.use_zhou_bhatt()) min_visit *= 2;
      histogram.clear();
      mcs.next_stage();
    } else {
      std::cerr << "stage " << mcs.stage() << ": histogram is not flat yet\n";
      histogram.subtract();
      mcs.reset_stage();
    }
  }
  if (mcs.doing_multicanonical() && mcs() == mcs.mcs_sweeps())
    histogram.store(measurements, "Partition Function Coefficient", "Histogram",
                    mcs.doing_multicanonical());
}


//
// diagonal update and cluster construction
//

void loop_worker::build()
{
  // initialize spin & operator information
  int nop = operators.size();
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(vlattice());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double bw = total_graph_weight();
  bool try_gap = true;
  for (operator_iterator opi = operators_p.begin();
       try_gap || opi != operators_p.end();) {

    // diagonal update & labeling
    if (try_gap) {
      if ((nop+1) * random() < bw * histogram.accept_rate(nop, nop+1)) {
        loop_graph_t g = choose_graph();
        if ((is_bond(g) &&
             is_compatible(g, spins_c[vsource(pos(g), vlattice())],
                           spins_c[vtarget(pos(g), vlattice())])) ||
            (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
          operators.push_back(local_operator_t(g));
          ++nop;
        } else {
          try_gap = false;
          histogram.visit(nop, logf, !mcs.doing_multicanonical());
          continue;
        }
      } else {
        try_gap = false;
        histogram.visit(nop, logf, !mcs.doing_multicanonical());
        continue;
      }
    } else {
      if (opi->is_diagonal()) {
        if (bw * random() < nop * histogram.accept_rate(nop, nop-1)) {
          --nop;
          ++opi;
          histogram.visit(nop, logf, !mcs.doing_multicanonical());
          continue;
        } else {
          if (opi->is_site()) {
            opi->assign_graph(choose_diagonal(opi->loc(),
              spins_c[opi->pos()]));
          } else {
            opi->assign_graph(choose_diagonal(opi->loc(),
              spins_c[vsource(opi->pos(), vlattice())],
              spins_c[vtarget(opi->pos(), vlattice())]));
          }
        }
      } else {
        opi->assign_graph(choose_offdiagonal(opi->loc()));
      }
      operators.push_back(*opi);
      ++opi;
      try_gap = true;
    }
    histogram.visit(nop, logf, !mcs.doing_multicanonical());

    operator_iterator oi = operators.end() - 1;
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice());
      int s1 = vtarget(oi->pos(), vlattice());
      boost::tie(current[s0], current[s1], oi->loop0, oi->loop1) =
        reconnect(fragments, oi->graph(), current[s0], current[s1]);
      if (oi->is_offdiagonal()) {
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
    } else {
      int s = oi->pos();
      boost::tie(current[s], oi->loop0, oi->loop1) =
        reconnect(fragments, oi->graph(), current[s]);
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
    }
  }

  // connect bottom and top cluster fragments after random permutation
  alps::fixed_capacity_vector<int, loop_config::max_2s> r;
  site_iterator rsi, rsi_end;
  for (boost::tie(rsi, rsi_end) = sites(rgraph()); rsi != rsi_end; ++rsi) {
    site_iterator vsi, vsi_end;
    boost::tie(vsi, vsi_end) = virtual_sites(vlattice(), rgraph(), *rsi);
    int offset = *vsi;
    int s2 = *vsi_end - *vsi;
    if (s2 == 1) {
      unify(fragments, offset, current[offset]);
    } else if (s2 > 1) {
      r.resize(s2);
      for (int i = 0; i < s2; ++i) r[i] = i;
      looper::restricted_random_shuffle(r.begin(), r.end(),
        spins.begin() + offset, spins_c.begin() + offset, random);
      for (int i = 0; i < s2; ++i)
        unify(fragments, offset+i, current[offset+r[i]]);
    }
  }
}


//
// cluster flip
//

template<typename BIPARTITE, typename FIELD, typename SIGN, typename IMPROVE>
void loop_worker::flip()
{
  if (!(is_bipartite() == BIPARTITE() &&
        has_field() == FIELD() &&
        is_signed() == SIGN() &&
        use_improved_estimator() == IMPROVE())) return;

  int nvs = num_sites(vlattice());
  int nop = operators.size();

  // assign cluster id
  int nc = 0;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) if (fi->is_root()) fi->id = nc++;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) fi->id = cluster_id(fragments, *fi);
  clusters.resize(0); clusters.resize(nc);

  std::copy(spins.begin(), spins.end(), spins_c.begin());
  if (IMPROVE()) estimates.resize(0); estimates.resize(nc);
  cluster_info_t::accumulator<cluster_fragment_t, FIELD, SIGN, IMPROVE>
    weight(clusters, fragments, field(), bond_sign(), site_sign());
  typename looper::measurement::accumulator<estimator_t, lattice_graph_t,
    time_t, cluster_fragment_t, BIPARTITE, IMPROVE>::type
    accum(estimates, fragments, vgraph());
  double t = 0;
  for (std::vector<local_operator_t>::iterator oi = operators.begin();
       oi != operators.end(); ++oi, t += 1) {
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice());
      int s1 = vtarget(oi->pos(), vlattice());
      weight.bond_sign(oi->loop_0(), oi->loop_1(), oi->pos());
      accum.term(oi->loop_l0(), t, s0, spins_c[s0]);
      accum.term(oi->loop_l1(), t, s1, spins_c[s1]);
      if (oi->is_offdiagonal()) {
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
      accum.start(oi->loop_u0(), t, s0, spins_c[s0]);
      accum.start(oi->loop_u1(), t, s1, spins_c[s1]);
    } else {
      int s = oi->pos();
      weight.site_sign(oi->loop_0(), oi->loop_1(), oi->pos());
      accum.term(oi->loop_l(), t, s, spins_c[s]);
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
      accum.start(oi->loop_u(), t, s, spins_c[s]);
    }
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    accum.start(s, 0, s, spins[s]);
    accum.term(current[s], nop, s, spins_c[s]);
    accum.at_zero(s, s, spins[s]);
  }

  // determine whether clusters are flipped or not
  double improved_sign = 1;
  for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
       ci != clusters.end(); ++ci) {
    ci->to_flip = (2*random()-1 < 0);
    if (SIGN() && IMPROVE()) if (ci->sign & 1 == 1) improved_sign = 0;
  }

  // flip operators & spins
  for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
    if (clusters[fragments[oi->loop_0()].id].to_flip ^
        clusters[fragments[oi->loop_1()].id].to_flip) oi->flip();
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id].to_flip) spins[s] ^= 1;

  // improved measurement
  if (IMPROVE()) {
    typename looper::measurement::collector<estimator_t, qmc_type, BIPARTITE,
      IMPROVE>::type coll;
    coll = std::accumulate(estimates.begin(), estimates.end(), coll);
    measurement_histograms.set_position(nop);
    coll.commit(measurement_histograms, beta(), num_sites(rgraph()), nop,
                improved_sign);
    if (SIGN()) measurement_histograms["Sign"] << improved_sign;
  }
}


//
// measurements
//

template<typename BIPARTITE, typename IMPROVE>
void loop_worker::measure()
{
  if (!(is_bipartite() == BIPARTITE() &&
        use_improved_estimator() == IMPROVE())) return;
  if (!mcs.doing_multicanonical()) return;

  int nrs = num_sites(rgraph());

  // sign
  double sign = 1;
  if (is_signed()) {
    int n = 0;
    for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
      if (oi->is_offdiagonal())
        n += (oi->is_bond()) ? bond_sign(oi->pos()) : site_sign(oi->pos());
    if (n & 1 == 1) sign = -1;
    if (!use_improved_estimator()) measurements["Sign"] << sign;
  }

  // energy
  int nop = operators.size();

  measurement_histograms.set_position(nop);
  looper::measurement::normal_estimator<estimator_t, qmc_type, BIPARTITE,
    IMPROVE>::type::measure(measurements, vgraph(), beta(), nrs, nop, sign,
                            spins, operators, spins_c);
}


//
// dynamic registration to the factories
//

const bool worker_registered =
  loop_factory::instance()->register_worker<loop_worker>("SSE QWL");

} // end namespace
