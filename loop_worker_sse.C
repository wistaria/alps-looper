/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_worker_sse.h"
#include <looper/permutation.h>
#include <alps/fixed_capacity_vector.h>
#include <boost/mpl/bool.hpp>

qmc_worker_sse::qmc_worker_sse(alps::ProcessList const& w,
                               alps::Parameters const& p, int n)
  : super_type(w, p, n, looper::is_path_integral<qmc_type>::type())
{
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

  using alps::RealObservable;
  using alps::make_observable;

  if (is_signed()) measurements << RealObservable("Sign");

  measurements
    << make_observable(
         RealObservable("Energy"), is_signed())
    << make_observable(
         RealObservable("Energy Density"), is_signed())
    << make_observable(
         RealObservable("Energy^2"), is_signed());

  normal_estimator_t::init(is_bipartite(), is_signed(), measurements);
  if (use_improved_estimator())
    improved_estimator_t::init(is_bipartite(), is_signed(), measurements);
}

void qmc_worker_sse::dostep()
{
  if (!can_work()) return;
  super_type::dostep();

  namespace mpl = boost::mpl;
  //          BIPARTITE    FIELD        SIGN         IMPROVE
  dostep_impl<mpl::true_,  mpl::true_,  mpl::true_,  mpl::true_ >();
  dostep_impl<mpl::true_,  mpl::true_,  mpl::true_,  mpl::false_>();
  dostep_impl<mpl::true_,  mpl::true_,  mpl::false_, mpl::true_ >();
  dostep_impl<mpl::true_,  mpl::true_,  mpl::false_, mpl::false_>();
  dostep_impl<mpl::true_,  mpl::false_, mpl::true_,  mpl::true_ >();
  dostep_impl<mpl::true_,  mpl::false_, mpl::true_,  mpl::false_>();
  dostep_impl<mpl::true_,  mpl::false_, mpl::false_, mpl::true_ >();
  dostep_impl<mpl::true_,  mpl::false_, mpl::false_, mpl::false_>();
  dostep_impl<mpl::false_, mpl::true_,  mpl::true_,  mpl::true_ >();
  dostep_impl<mpl::false_, mpl::true_,  mpl::true_,  mpl::true_ >();
  dostep_impl<mpl::false_, mpl::true_,  mpl::false_, mpl::true_ >();
  dostep_impl<mpl::false_, mpl::true_,  mpl::false_, mpl::false_>();
  dostep_impl<mpl::false_, mpl::false_, mpl::true_,  mpl::true_ >();
  dostep_impl<mpl::false_, mpl::false_, mpl::true_,  mpl::true_ >();
  dostep_impl<mpl::false_, mpl::false_, mpl::false_, mpl::true_ >();
  dostep_impl<mpl::false_, mpl::false_, mpl::false_, mpl::false_>();
}

template<typename BIPARTITE, typename FIELD, typename SIGN, typename IMPROVE>
void qmc_worker_sse::dostep_impl()
{
  if (!(is_bipartite() == BIPARTITE() &&
        has_field() == FIELD() &&
        is_signed() == SIGN() &&
        use_improved_estimator() == IMPROVE())) return;

  //
  // diagonal update and cluster construction
  //

  // initialize spin & operator information
  int nop = operators.size();
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(vgraph());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;

  double bw = beta() * total_graph_weight();
  bool try_gap = true;
  for (operator_iterator opi = operators_p.begin();
       try_gap || opi != operators_p.end();) {

    // diagonal update & labeling
    if (try_gap) {
      if ((nop + 1) * random() < bw) {
        loop_graph_t g = choose_graph();
        if ((is_bond(g) &&
             is_compatible(g, spins_c[vsource(pos(g), vlattice())],
                           spins_c[vtarget(pos(g), vlattice())])) ||
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

  //
  // cluster flip
  //

  // assign cluster id
  int nc = 0;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) if (fi->is_root()) fi->id = nc++;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) fi->id = cluster_id(fragments, *fi);
  clusters.resize(0); clusters.resize(nc);

  std::copy(spins.begin(), spins.end(), spins_c.begin());
  if (IMPROVE()) estimates.resize(0); estimates.resize(nc);
  cluster_info_t::accumulator<cluster_fragment_t, FIELD>
    weight(clusters, fragments, field());
  typename improved_estimator_t::accumulator<lattice_graph_t,
    cluster_fragment_t, IMPROVE, BIPARTITE>::type
    accum(estimates, fragments, vgraph());
  double t = 0;
  for (std::vector<local_operator_t>::iterator oi = operators.begin();
       oi != operators.end(); ++oi, t += 1) {
    if (oi->is_bond()) {
      int s0 = vsource(oi->pos(), vlattice());
      int s1 = vtarget(oi->pos(), vlattice());
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
      weight.term(oi->loop_l(), t, s, spins_c[s]);
      accum.term(oi->loop_l(), t, s, spins_c[s]);
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
      weight.start(oi->loop_u(), t, s, spins_c[s]);
      accum.start(oi->loop_u(), t, s, spins_c[s]);
    }
  }
  int nop_or_1 = std::max(nop, 1);
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.start(s, 0, s, spins[s]);
    weight.term(current[s], nop_or_1, s, spins_c[s]);
    accum.start(s, 0, s, spins[s]);
    accum.term(current[s], nop, s, spins_c[s]);
    accum.at_zero(s, s, spins[s]);
  }

  // determine whether clusters are flipped or not
  for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
       ci != clusters.end(); ++ci)
    ci->to_flip =
      ((2*random()-1) < (FIELD() ? std::tanh(beta()*ci->weight/nop_or_1) : 0));

  // flip operators & spins
  for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
    if (clusters[fragments[oi->loop_0()].id].to_flip ^
        clusters[fragments[oi->loop_1()].id].to_flip) oi->flip();
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id].to_flip) spins[s] ^= 1;

  //
  // measurements
  //

  int nrs = num_sites(rgraph());

  // sign
  double sign = 1;
  if (SIGN()) {
    int n = 0;
    for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
      if (oi->is_offdiagonal())
        if (oi->is_bond())
          n ^= bond_sign(oi->pos());
        else
          n ^= site_sign(oi->pos());
    sign = 1 - 2 * n;
    /* if (!IMPROVE()) */ measurements["Sign"] << sign;
  }

  // energy
  double ene = energy_offset() - nop / beta();
  measurements["Energy"] << sign * ene;
  measurements["Energy Density"] << sign * ene / nrs;
  measurements["Energy^2"]
    << sign * (looper::power2(ene) - nop / looper::power2(beta()));

  // other measurements
  if (IMPROVE()) {
    std::accumulate(estimates.begin(), estimates.end(),
      typename improved_estimator_t::collector<BIPARTITE>::type()).
      commit(qmc_type(), nrs, beta(), nop, measurements);
  } else {
    normal_estimator_t::do_measurement(qmc_type(), vgraph(), BIPARTITE(), nrs,
      beta(), nop, sign, spins, operators, spins_c, measurements);
  }
}

void qmc_worker_sse::save(alps::ODump& dp) const
{
  super_type::save(dp);
  dp << spins << operators;
}

void qmc_worker_sse::load(alps::IDump& dp)
{
  super_type::load(dp);
  dp >> spins >> operators;
}
