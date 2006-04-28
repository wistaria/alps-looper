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

#include "loop_config.h"
#include "loop_factory.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/model.h>
#include <looper/montecarlo.h>
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/type.h>
#include <looper/weight.h>
#include <alps/fixed_capacity_vector.h>

namespace {

class loop_worker
  : public alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
{
public:
  typedef looper::path_integral                            qmc_type;
  typedef alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
                                                           super_type;

  typedef loop_config::lattice_graph_t                     lattice_graph_t;
  typedef loop_config::time_t                              time_t;
  typedef loop_config::loop_graph_t                        loop_graph_t;
  typedef looper::local_operator<qmc_type, loop_graph_t, time_t>
                                                           local_operator_t;
  typedef std::vector<local_operator_t>                    operator_string_t;
  typedef operator_string_t::iterator                      operator_iterator;

  typedef looper::union_find::node                         cluster_fragment_t;
  typedef looper::cluster_info                             cluster_info_t;

  typedef loop_config::estimator_t                         estimator_t;
  typedef looper::measurement::estimate<estimator_t>::type estimate_t;

  loop_worker(alps::ProcessList const& w, alps::Parameters const& p, int n);
  void dostep();

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double work_done() const { return mcs.progress(); }

  void save(alps::ODump& dp) const
  { super_type::save(dp); dp << mcs << spins << operators; }
  void load(alps::IDump& dp)
  { super_type::load(dp); dp >> mcs >> spins >> operators; }

protected:
  void build();
  template<typename BIPARTITE, typename FIELD, typename SIGN, typename IMPROVE>
  void flip();
  template<typename BIPARTITE, typename IMPROVE>
  void measure();

private:
  // parameters
  double beta;
  double energy_offset;
  bool has_field, is_frustrated, is_signed, use_improved_estimator;
  std::vector<double> field;
  std::vector<int> bond_sign, site_sign;
  looper::virtual_lattice<lattice_graph_t> vlattice;

  // random number generator
  looper::graph_chooser<loop_graph_t, super_type::engine_type> chooser;
  boost::variate_generator<super_type::engine_type&,
    boost::exponential_distribution<> > r_time;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // working vectors
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
  : super_type(w, p, n),
    chooser(*engine_ptr),
    r_time(*engine_ptr, boost::exponential_distribution<>()),
    mcs(p)
{
  beta = 1.0 / alps::evaluate("T", p);
  if (beta < 0)
    boost::throw_exception(std::invalid_argument("negative temperature"));

  looper::model_parameter mp(p, *this);
  energy_offset = mp.energy_offset();
  has_field = mp.has_field();
  is_frustrated = mp.is_frustrated();
  is_signed = mp.is_signed();
  use_improved_estimator =
    !has_field && !p.defined("DISABLE_IMPROVED_ESTIMATOR");

  if (has_field) std::cerr << "WARNING: model has magnetic field\n";
  if (is_frustrated)
    std::cerr << "WARNING: model is classically frustrated\n";
  if (is_signed) std::cerr << "WARNING: model has negative signs\n";
  if (!use_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  vlattice.generate(graph(), mp, mp.has_d_term());

  double fs = p.value_or_default("FORCE_SCATTER", is_frustrated ? 0.1 : 0);
  looper::weight_table wt(mp, graph(), vlattice, fs);
  energy_offset += wt.energy_offset();
  chooser.init(wt, looper::is_path_integral<qmc_type>::type());
  r_time.distribution() =
    boost::exponential_distribution<>(beta * chooser.weight());

  if (is_signed) {
    bond_sign.resize(num_bonds(vlattice));
    looper::weight_table::bond_weight_iterator bi, bi_end;
    for (boost::tie(bi, bi_end) = wt.bond_weights(); bi != bi_end; ++bi)
      bond_sign[bi->first] = (bi->second.sign < 0) ? 1 : 0;
    site_sign.resize(num_sites(vlattice));
    looper::weight_table::site_weight_iterator si, si_end;
    for (boost::tie(si, si_end) = wt.site_weights(); si != si_end; ++si)
      site_sign[si->first] = (si->second.sign < 0) ? 1 : 0;
  }

  if (has_field) {
    field.resize(num_sites(vlattice));
    int i = 0;
    site_iterator rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(); rsi != rsi_end; ++rsi) {
      site_iterator vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) = virtual_sites(vlattice, graph(), *rsi);
           vsi != vsi_end; ++vsi, ++i)
        field[i] = mp.site(*rsi, graph()).hz;
    }
  }

  // configuration
  int nvs = num_sites(vlattice);
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);

  // measurements
  measurements <<
    make_observable(alps::SimpleRealObservable("Temperature"));
  measurements <<
    make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  measurements <<
    make_observable(alps::SimpleRealObservable("Number of Sites"));
  measurements <<
    make_observable(alps::SimpleRealObservable("Number of Clusters"));
  if (is_signed) measurements << alps::RealObservable("Sign");
  looper::energy_estimator::initialize(measurements, is_signed);
  estimator_t::initialize(measurements, is_bipartite(), is_signed,
                          use_improved_estimator);
}

void loop_worker::dostep()
{
  namespace mpl = boost::mpl;

  if (!mcs.can_work()) return;
  ++mcs;

  build();

  //   BIPARTITE    FIELD        SIGN         IMPROVE
  flip<mpl::true_,  mpl::true_,  mpl::true_,  mpl::true_ >();
  flip<mpl::true_,  mpl::true_,  mpl::true_,  mpl::false_>();
  flip<mpl::true_,  mpl::true_,  mpl::false_, mpl::true_ >();
  flip<mpl::true_,  mpl::true_,  mpl::false_, mpl::false_>();
  flip<mpl::true_,  mpl::false_, mpl::true_,  mpl::true_ >();
  flip<mpl::true_,  mpl::false_, mpl::true_,  mpl::false_>();
  flip<mpl::true_,  mpl::false_, mpl::false_, mpl::true_ >();
  flip<mpl::true_,  mpl::false_, mpl::false_, mpl::false_>();
  flip<mpl::false_, mpl::true_,  mpl::true_,  mpl::true_ >();
  flip<mpl::false_, mpl::true_,  mpl::true_,  mpl::false_>();
  flip<mpl::false_, mpl::true_,  mpl::false_, mpl::true_ >();
  flip<mpl::false_, mpl::true_,  mpl::false_, mpl::false_>();
  flip<mpl::false_, mpl::false_, mpl::true_,  mpl::true_ >();
  flip<mpl::false_, mpl::false_, mpl::true_,  mpl::false_>();
  flip<mpl::false_, mpl::false_, mpl::false_, mpl::true_ >();
  flip<mpl::false_, mpl::false_, mpl::false_, mpl::false_>();

  //      BIPARTITE    IMPROVE
  measure<mpl::true_,  mpl::true_ >();
  measure<mpl::true_,  mpl::false_>();
  measure<mpl::false_, mpl::true_ >();
  measure<mpl::false_, mpl::false_>();
}


//
// diagonal update and cluster construction
//

void loop_worker::build()
{
  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(vlattice);
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double t = r_time();
  for (operator_iterator opi = operators_p.begin();
       t < 1 || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      loop_graph_t g = chooser.graph();
      if (((is_bond(g) &&
            is_compatible(g, spins_c[vsource(pos(g), vlattice)],
                          spins_c[vtarget(pos(g), vlattice)])) ||
          (is_site(g) && is_compatible(g, spins_c[pos(g)])))) {
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
        opi->assign_graph(chooser.offdiagonal(opi->loc()));
        operators.push_back(*opi);
        ++opi;
      }
    }

    operator_iterator oi = operators.end() - 1;
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice);
      int s1 = vtarget(oi->pos(), vlattice);
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

  // symmetrize spins
  alps::fixed_capacity_vector<int, loop_config::max_2s> r;
  site_iterator rsi, rsi_end;
  for (boost::tie(rsi, rsi_end) = sites(); rsi != rsi_end; ++rsi) {
    site_iterator vsi, vsi_end;
    boost::tie(vsi, vsi_end) = virtual_sites(vlattice, graph(), *rsi);
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
        has_field == FIELD() &&
        is_signed == SIGN() &&
        use_improved_estimator == IMPROVE())) return;

  int nvs = num_sites(vlattice);
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
    weight(clusters, fragments, field, bond_sign, site_sign);
  typename looper::measurement::accumulator<estimator_t, lattice_graph_t,
    time_t, cluster_fragment_t, BIPARTITE, IMPROVE>::type
    accum(estimates, fragments, vlattice.graph());
  for (std::vector<local_operator_t>::iterator oi = operators.begin();
       oi != operators.end(); ++oi) {
    time_t t = oi->time();
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice);
      int s1 = vtarget(oi->pos(), vlattice);
      weight.bond_sign(oi->loop_0(), oi->loop_1(), oi->pos());
      weight.term(oi->loop_l0(), t, s0, spins_c[s0]);
      weight.term(oi->loop_l1(), t, s1, spins_c[s1]);
      accum.term(oi->loop_l0(), t, s0, spins_c[s0]);
      accum.term(oi->loop_l1(), t, s1, spins_c[s1]);
      if (oi->is_offdiagonal()) {
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
      weight.start(oi->loop_u0(), t, s0, spins_c[s0]);
      weight.start(oi->loop_u1(), t, s1, spins_c[s1]);
      accum.start(oi->loop_u0(), t, s0, spins_c[s0]);
      accum.start(oi->loop_u1(), t, s1, spins_c[s1]);
    } else {
      int s = oi->pos();
      weight.site_sign(oi->loop_0(), oi->loop_1(), oi->pos());
      weight.term(oi->loop_l(), t, s, spins_c[s]);
      accum.term(oi->loop_l(), t, s, spins_c[s]);
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
      weight.start(oi->loop_u(), t, s, spins_c[s]);
      accum.start(oi->loop_u(), t, s, spins_c[s]);
    }
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.start(s, time_t(0), s, spins[s]);
    weight.term(current[s], time_t(1), s, spins_c[s]);
    accum.start(s, time_t(0), s, spins[s]);
    accum.term(current[s], time_t(1), s, spins_c[s]);
    accum.at_zero(s, s, spins[s]);
  }

  // determine whether clusters are flipped or not
  double improved_sign = 1;
  for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
       ci != clusters.end(); ++ci) {
    ci->to_flip = ((2*random()-1) < (FIELD() ? std::tanh(ci->weight) : 0));
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
    coll.commit(measurements, beta, num_sites(), nop, improved_sign);
    if (SIGN()) measurements["Sign"] << improved_sign;
  }
  measurements["Number of Clusters"] << (double)clusters.size();
}


//
// measurements
//

template<typename BIPARTITE, typename IMPROVE>
void loop_worker::measure()
{
  if (!(is_bipartite() == BIPARTITE() &&
        use_improved_estimator == IMPROVE())) return;

  int nrs = num_sites();
  measurements["Temperature"] << 1/beta;
  measurements["Inverse Temperature"] << beta;
  measurements["Number of Sites"] << (double)nrs;

  // sign
  double sign = 1;
  if (is_signed) {
    int n = 0;
    for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
      if (oi->is_offdiagonal())
        n += (oi->is_bond()) ? bond_sign[oi->pos()] : site_sign[oi->pos()];
    if (n & 1 == 1) sign = -1;
    if (!use_improved_estimator) measurements["Sign"] << sign;
  }

  // energy
  int nop = operators.size();
  double ene = energy_offset - nop / beta;
  if (has_field) {
    for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
         ci != clusters.end(); ++ci)
      if (ci->to_flip) ene -= ci->weight;
      else ene += ci->weight;
  }
  looper::energy_estimator::measure(measurements, beta, nrs, nop, sign, ene);

  // other quantities
  looper::measurement::normal_estimator<estimator_t, qmc_type, BIPARTITE,
    IMPROVE>::type::measure(measurements, vlattice.graph(), beta, nrs, nop,
                            sign, spins, operators, spins_c);
}


//
// dynamic registration to the factories
//

const bool worker_registered =
  loop_factory::instance()->register_worker<loop_worker>("path integral");
const bool evaluator_registered = evaluator_factory::instance()->
  register_evaluator<looper::evaluator<loop_config::estimator_t> >
  ("path integral");

} // end namespace
