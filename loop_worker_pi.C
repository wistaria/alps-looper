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
#include <looper/operator.h>
#include <looper/permutation.h>
#include <looper/type.h>
#include <alps/fixed_capacity_vector.h>

class qmc_worker_pi : public qmc_worker
{
public:
  typedef looper::path_integral         qmc_type;
  typedef qmc_worker                    super_type;

  typedef looper::local_operator<qmc_type, loop_graph_t, time_t>
                                        local_operator_t;
  typedef std::vector<local_operator_t> operator_string_t;
  typedef operator_string_t::iterator   operator_iterator;

  typedef looper::union_find::node      cluster_fragment_t;
  typedef looper::cluster_info          cluster_info_t;

  typedef loop_config::estimator_t      estimator_t;
  typedef looper::measurement::estimate<estimator_t>::type estimate_t;

  qmc_worker_pi(alps::ProcessList const& w, alps::Parameters const& p, int n);
  virtual void dostep();

  void save(alps::ODump& dp) const
  { super_type::save(dp); dp << spins << operators; }
  void load(alps::IDump& dp)
  { super_type::load(dp); dp >> spins >> operators; }

protected:
  void build();
  template<typename BIPARTITE, typename FIELD, typename SIGN, typename IMPROVE>
  void flip();
  template<typename BIPARTITE, typename IMPROVE>
  void measure();

private:
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // working vectors (no checkpointing)
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<cluster_info_t> clusters;
  std::vector<estimate_t> estimates;
};

qmc_worker_pi::qmc_worker_pi(alps::ProcessList const& w,
                             alps::Parameters const& p, int n)
  : super_type(w, p, n, looper::is_path_integral<qmc_type>::type())
{
  if (w == alps::ProcessList()) return;

  //
  // initialize configuration
  //

  int nvs = num_sites(vgraph());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
  operators.resize(0);
  spins_c.resize(nvs);
  current.resize(nvs);

  //
  // init measurements
  //

  if (is_signed()) measurements << alps::RealObservable("Sign");
  estimator_t::initialize(measurements, is_bipartite(), is_signed(),
                          use_improved_estimator());
}

void qmc_worker_pi::dostep()
{
  namespace mpl = boost::mpl;

  if (!can_work()) return;
  super_type::dostep();

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
  flip<mpl::false_, mpl::true_,  mpl::true_,  mpl::true_ >();
  flip<mpl::false_, mpl::true_,  mpl::false_, mpl::true_ >();
  flip<mpl::false_, mpl::true_,  mpl::false_, mpl::false_>();
  flip<mpl::false_, mpl::false_, mpl::true_,  mpl::true_ >();
  flip<mpl::false_, mpl::false_, mpl::true_,  mpl::true_ >();
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

void qmc_worker_pi::build()
{
  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(vgraph());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double t = advance();
  for (operator_iterator opi = operators_p.begin();
       t < 1 || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      loop_graph_t g = choose_graph();
      if (((is_bond(g) &&
            is_compatible(g, spins_c[vsource(pos(g), vlattice())],
                          spins_c[vtarget(pos(g), vlattice())])) ||
          (is_site(g) && is_compatible(g, spins_c[pos(g)])))) {
        operators.push_back(local_operator_t(g, t));
        t += advance();
      } else {
        t += advance();
        continue;
      }
    } else {
      if (opi->is_diagonal()) {
        ++opi;
        continue;
      } else {
        opi->assign_graph(choose_offdiagonal(opi->loc()));
        operators.push_back(*opi);
        ++opi;
      }
    }

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
void qmc_worker_pi::flip()
{
  if (!(is_bipartite() == BIPARTITE() &&
        has_field() == FIELD() &&
        is_signed() == SIGN() &&
        use_improved_estimator() == IMPROVE())) return;

  int nvs = num_sites(vgraph());
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
  for (std::vector<local_operator_t>::iterator oi = operators.begin();
       oi != operators.end(); ++oi) {
    time_t t = oi->time();
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice());
      int s1 = vtarget(oi->pos(), vlattice());
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
    coll.commit(measurements, beta(), num_sites(rgraph()), nop, improved_sign);
    if (SIGN()) measurements["Sign"] << improved_sign;
  }
}


//
// measurements
//

template<typename BIPARTITE, typename IMPROVE>
void qmc_worker_pi::measure()
{
  if (!(is_bipartite() == BIPARTITE() &&
        use_improved_estimator() == IMPROVE())) return;

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
  double ene = energy_offset() - nop / beta();
  if (has_field()) {
    for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
         ci != clusters.end(); ++ci)
      if (ci->to_flip) ene -= ci->weight;
      else ene += ci->weight;
  }

  looper::measurement::normal_estimator<estimator_t, qmc_type, BIPARTITE,
    IMPROVE>::type::measure(measurements, vgraph(), beta(), nrs, nop, sign,
                            ene, spins, operators, spins_c);
}

//
// dynamic registration to the qmc_factory
//

namespace {

const bool registered =
  qmc_factory::instance()->register_worker<qmc_worker_pi>("path integral");

}
