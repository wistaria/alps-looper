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

#include "loop_worker_pi.h"
#include <looper/permutation.h>
#include <alps/fixed_capacity_vector.h>

qmc_worker_pi::qmc_worker_pi(const alps::ProcessList& w,
                             const alps::Parameters& p, int n)
  : super_type(w, p, n, true)
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

  measurements
    << make_observable(
         RealObservable("Magnetization"), is_signed())
    << make_observable(
         RealObservable("Magnetization^2"), is_signed())
    << make_observable(
         RealObservable("Susceptibility"), is_signed());

  if (is_bipartite()) {
    measurements
      << make_observable(
           RealObservable("Staggered Magnetization"), is_signed())
      << make_observable(
           RealObservable("Staggered Magnetization^2"), is_signed())
      << make_observable(
           RealObservable("Staggered Susceptibility"), is_signed());
  }

  if (!is_classically_frustrated()) {
    measurements
      << RealObservable("Generalized Magnetization^2")
      << RealObservable("Generalized Susceptibility");
    if (is_bipartite()) {
      measurements
        << RealObservable("Staggered Generalized Magnetization^2")
        << RealObservable("Staggered Generalized Susceptibility");
    }
  }
}

void qmc_worker_pi::dostep()
{
  if (!can_work()) return;
  super_type::dostep();

  //
  // diagonal update and cluster construction
  //

  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_c.begin());
  std::swap(operators, operators_p); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(vgraph());
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;
  int ghost = has_longitudinal_field() ? add(fragments) : 0;

  double t = advance();
  for (operator_iterator opi = operators_p.begin();
       t < beta() || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      // insert diagonal operator and graph if compatible
      loop_graph_t g = choose_graph();
      if (
          ((is_bond(g) &&
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
        // remove diagonal operator with probability one (= nothing to do)
        ++opi;
        continue;
      } else {
        // assign graph to offdiagonal operator
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
      if (oi->is_locked()) unify(fragments, ghost, current[s]);
      if (oi->is_offdiagonal()) spins_c[s] ^= 1;
    }
  }

  // connect bottom and top cluster fragments after random permutation
  {
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

  // assign cluster id & determine if clusters are to be flipped
  clusters.resize(0);
  for (std::vector<cluster_fragment_t>::iterator ci = fragments.begin();
       ci != fragments.end(); ++ci)
    if (ci->is_root()) {
      ci->id = clusters.size();
      clusters.push_back(cluster_info_t(random() < 0.5));
    }
  if (has_longitudinal_field())
    clusters[cluster_id(fragments, ghost)].to_flip = false;

  // flip operators
  for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
    if (clusters[cluster_id(fragments, oi->loop0)].to_flip ^
        clusters[cluster_id(fragments, oi->loop1)].to_flip) oi->flip();
  for (int s = 0; s < nvs; ++s)
    if (clusters[cluster_id(fragments, s)].to_flip) spins[s] ^= 1;

  //
  // measurements
  //

  using std::sqrt;
  using looper::sqr;
  double nrsi = 1.0 / (double)num_sites(rgraph());

  // energy
  int nop = operators.size();
  double ene = energy_offset() - nop / beta();
  measurements["Energy"] << ene;
  measurements["Energy Density"] << nrsi * ene;
  measurements["Energy^2"] << sqr(ene) - nop / sqr(beta());
  // measurements["beta * Energy / sqrt(N)"] << sqrt(nrsi) * beta * ene;
  // measurements["beta * Energy^2"] <<
  //   nrsi * sqr(beta) * (sqr(ene) - nop / sqr(beta));

  // magnetization && susceptibility
  if (is_bipartite()) {
    int nm = 0;
    int ns = 0;
    site_iterator si, si_end;
    for (boost::tie(si, si_end) = sites(vgraph()); si != si_end; ++si) {
      nm += (1 - 2 * spins[*si]);
      ns += looper::gauge(vgraph(), *si) * (1 - 2 * spins[*si]);
    }
    double umag = nm / 2.0;
    double smag = ns / 2.0;
    measurements["Magnetization"] << nrsi * umag;
    measurements["Magnetization^2"] << nrsi * umag * umag;
    measurements["Staggered Magnetization"] << nrsi * smag;
    measurements["Staggered Magnetization^2"] << nrsi * smag * smag;

    double umag_a = 0; /* 0 * umag; */
    double smag_a = 0; /* 0 * smag; */
    std::copy(spins.begin(), spins.end(), spins_c.begin());
    for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
      if (oi->is_offdiagonal()) {
         umag_a += oi->time() * umag;
         smag_a += oi->time() * smag;
         if (oi->is_site()) {
           unsigned int s = oi->pos();
           spins_c[s] ^= 1;
           umag += (1 - 2 * spins_c[s]);
          smag += looper::gauge(vgraph(), s) * (1 - 2 * spins_c[s]);
         } else {
           unsigned int s0 = vsource(oi->pos(), vlattice());
           unsigned int s1 = vtarget(oi->pos(), vlattice());
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
           umag += 1 - 2 * spins_c[s0] + 1 - 2 * spins_c[s1];
          smag += looper::gauge(vgraph(), s0) * (1 - 2 * spins_c[s0])
            + looper::gauge(vgraph(), s1) * (1 - 2 * spins_c[s1]);
         }
         umag_a -= oi->time() * umag;
         smag_a -= oi->time() * smag;
      }
    umag_a += beta() * umag;
    smag_a += beta() * smag;
    measurements["Susceptibility"]
      << nrsi * umag_a * umag_a / beta();
    measurements["Staggered Susceptibility"]
      << nrsi * smag_a * smag_a / beta();
  } else {
    int nm = 0;
    site_iterator si, si_end;
    for (boost::tie(si, si_end) = sites(vgraph()); si != si_end; ++si) {
      nm += (1 - 2 * spins[*si]);
    }
    double umag = nm / 2.0;
    measurements["Magnetization"] << nrsi * umag;
    measurements["Magnetization^2"] << nrsi * umag * umag;

    double umag_a = 0 * umag;
    std::copy(spins.begin(), spins.end(), spins_c.begin());
    for (operator_iterator oi = operators.begin(); oi != operators.end(); ++oi)
      if (oi->is_offdiagonal()) {
        umag_a += oi->time() * umag;
        if (oi->is_site()) {
          unsigned int s = oi->pos();
          spins_c[s] ^= 1;
          umag += (1 - 2 * spins_c[s]);
        } else {
          unsigned int s0 = vsource(oi->pos(), vlattice());
          unsigned int s1 = vtarget(oi->pos(), vlattice());
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
          umag += 1 - 2 * spins_c[s0] + 1 - 2 * spins_c[s1];
        }
        umag_a -= oi->time() * umag;
      }
    umag_a += beta() * umag;
    measurements["Susceptibility"] << nrsi * umag_a * umag_a / beta();
  }
}

void qmc_worker_pi::save(alps::ODump& dp) const
{
  super_type::save(dp);
  dp << spins << operators;
}

void qmc_worker_pi::load(alps::IDump& dp)
{
  super_type::load(dp);
  dp >> spins >> operators;
}
