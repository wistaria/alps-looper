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

namespace {

class loop_worker
  : public alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
{
public:
  typedef looper::sse                                      qmc_type;
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
  bool is_frustrated, is_signed, use_improved_estimator;
  std::vector<double> field;
  std::vector<int> bond_sign, site_sign;
  looper::virtual_lattice<lattice_graph_t> vlattice;

  // random number generator
  looper::graph_chooser<loop_graph_t, super_type::engine_type> chooser;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // observables
  alps::ObservableSet& obs;

  // working vectors
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<cluster_info_t> clusters;
  std::vector<estimate_t> estimates;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::ProcessList const& w,
                         alps::Parameters const& p, int n)
  : super_type(w, p, n), chooser(*engine_ptr), mcs(p), obs(measurements)
{
  beta = 1.0 / alps::evaluate("T", p);
  if (beta < 0)
    boost::throw_exception(std::invalid_argument("negative temperature"));

  looper::model_parameter mp(p, *this);
  energy_offset = mp.energy_offset();
  is_frustrated = mp.is_frustrated();
  is_signed = mp.is_signed();
  use_improved_estimator = !p.defined("DISABLE_IMPROVED_ESTIMATOR");

  if (mp.has_field())
    boost::throw_exception(std::logic_error("longitudinal field is currently "
      "not supported in SSE representation"));
  if (is_frustrated)
    std::cerr << "WARNING: model is classically frustrated\n";
  if (is_signed) std::cerr << "WARNING: model has negative signs\n";
  if (!use_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  vlattice.generate(graph(), mp, mp.has_d_term());
  perm.resize(max_virtual_vertices(vlattice));

  double fs = p.value_or_default("FORCE_SCATTER", is_frustrated ? 0.1 : 0);
  looper::weight_table wt(mp, graph(), vlattice, fs);
  energy_offset += wt.energy_offset();
  chooser.init(wt, looper::is_path_integral<qmc_type>::type());

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

  // configuration
  int nvs = num_sites(vlattice);
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  spins_c.resize(nvs);
  current.resize(nvs);

  obs << make_observable(alps::SimpleRealObservable("Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Number of Sites"));
  obs << make_observable(alps::SimpleRealObservable("Number of Clusters"));
  if (is_signed) obs << alps::RealObservable("Sign");
  looper::energy_estimator::initialize(obs, is_signed);
  estimator_t::initialize(obs, is_bipartite(), is_signed,
                          use_improved_estimator);
}

void loop_worker::dostep()
{
  namespace mpl = boost::mpl;

  if (!mcs.can_work()) return;
  ++mcs;

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
  int nvs = num_sites(vlattice);
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double bw = beta * chooser.weight();
  bool try_gap = true;
  for (operator_iterator opi = operators_p.begin();
       try_gap || opi != operators_p.end();) {

    // diagonal update & labeling
    if (try_gap) {
      if ((nop+1) * random() < bw) {
        loop_graph_t g = chooser.graph();
        if ((is_bond(g) &&
             is_compatible(g, spins_c[vsource(pos(g), vlattice)],
                           spins_c[vtarget(pos(g), vlattice)])) ||
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
        if (bw * random() < nop) {
          --nop;
          ++opi;
          continue;
        } else {
          if (opi->is_site()) {
            opi->assign_graph(chooser.diagonal(opi->loc(),
              spins_c[opi->pos()]));
          } else {
            opi->assign_graph(chooser.diagonal(opi->loc(),
              spins_c[vsource(opi->pos(), vlattice)],
              spins_c[vtarget(opi->pos(), vlattice)]));
          }
        }
      } else {
        opi->assign_graph(chooser.offdiagonal(opi->loc()));
      }
      operators.push_back(*opi);
      ++opi;
      try_gap = true;
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
  if (max_virtual_vertices(vlattice) == 1) {
    for (int i = 0; i < nvs; ++i) unify(fragments, i, current[i]);
  } else {
    site_iterator rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(); rsi != rsi_end; ++rsi) {
      site_iterator vsi, vsi_end;
      boost::tie(vsi, vsi_end) = virtual_sites(vlattice, graph(), *rsi);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins_c.begin() + offset, random);
      for (int i = 0; i < s2; ++i)
        unify(fragments, offset+i, current[offset+perm[i]]);
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
  double t = 0;
  for (std::vector<local_operator_t>::iterator oi = operators.begin();
       oi != operators.end(); ++oi, t += 1) {
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice);
      int s1 = vtarget(oi->pos(), vlattice);
      weight.bond_sign(oi->loop_0(), oi->loop_1(), oi->pos());
      accum.term2(oi->loop_l0(), oi->loop_l1(), t, s0, s1,
                  spins_c[s0], spins_c[s1]);
      if (oi->is_offdiagonal()) {
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
      accum.start2(oi->loop_u0(), oi->loop_u1(), t, s0, s1,
                   spins_c[s0], spins_c[s1]);
    } else {
      int s = oi->pos();
      weight.site_sign(oi->loop_0(), oi->loop_1(), oi->pos());
      accum.term1(oi->loop_l(), t, s, spins_c[s]);
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
      accum.start1(oi->loop_u(), t, s, spins_c[s]);
    }
  }
  for (unsigned int s = 0; s < nvs; ++s) {
    accum.at_bot(s,          0,   s, spins[s]);
    accum.at_top(current[s], nop, s, spins_c[s]);
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
    coll.commit(obs, beta, num_sites(), nop, improved_sign);
    if (SIGN()) obs["Sign"] << improved_sign;
  }
  obs["Number of Clusters"] << (double)clusters.size();
}


//
// measurement
//

template<typename BIPARTITE, typename IMPROVE>
void loop_worker::measure()
{
  if (!(is_bipartite() == BIPARTITE() &&
        use_improved_estimator == IMPROVE())) return;

  int nrs = num_sites();
  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Number of Sites"] << (double)nrs;

  // sign
  double sign = 1;
  if (is_signed) {
    int n = 0;
    for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
      if (oi->is_offdiagonal())
        n += (oi->is_bond()) ? bond_sign[oi->pos()] : site_sign[oi->pos()];
    if (n & 1 == 1) sign = -1;
    if (!use_improved_estimator) obs["Sign"] << sign;
  }

  // energy
  int nop = operators.size();
  double ene = energy_offset - nop / beta;
  looper::energy_estimator::measure(obs, beta, nrs, nop, sign, ene);

  // other quantities
  looper::measurement::normal_estimator<estimator_t, qmc_type, BIPARTITE,
    IMPROVE>::type::measure(obs, vlattice.graph(), beta, nrs, nop, sign, spins,
                            operators, spins_c);
}


//
// dynamic registration to the factories
//

const bool worker_registered =
  loop_factory::instance()->register_worker<loop_worker>("SSE");
const bool evaluator_registered = evaluator_factory::instance()->
  register_evaluator<looper::evaluator<loop_config::estimator_t> >
  ("SSE");

} // end namespace
