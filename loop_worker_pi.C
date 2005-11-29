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
#include <alps/fixed_capacity_vector.h>

qmc_worker_pi::qmc_worker_pi(const alps::ProcessList& w,
                             const alps::Parameters& p, int n)
  : super_type(w, p, n),
    beta(1.0 / static_cast<double>(p["T"]))
{
  if (beta < 0)
    boost::throw_exception(std::invalid_argument("negative beta"));

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

  if (mp.is_signed()) {
    measurements << RealObservable("Sign");
  }

  measurements
    << make_observable(
         RealObservable("Energy"), mp.is_signed())
    << make_observable(
         RealObservable("Energy Density"), mp.is_signed())
    << make_observable(
         RealObservable("Diagonal Energy Density"), mp.is_signed())
    << make_observable(
         RealObservable("Energy Density^2"), mp.is_signed())
    << make_observable(
         RealObservable("beta * Energy / sqrt(N)"), mp.is_signed())
    << make_observable(
         RealObservable("beta * Energy^2"), mp.is_signed());

  measurements
    << make_observable(
         RealObservable("Magnetization"), mp.is_signed())
    << make_observable(
         RealObservable("Magnetization^2"), mp.is_signed())
    << make_observable(
         RealObservable("Susceptibility"), mp.is_signed());

  if (is_bipartite()) {
    measurements
      << make_observable(
           RealObservable("Staggered Magnetization"), mp.is_signed())
      << make_observable(
           RealObservable("Staggered Magnetization^2"), mp.is_signed())
      << make_observable(
           RealObservable("Staggered Susceptibility"), mp.is_signed());
  }

  if (!mp.is_classically_frustrated()) {
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
  int ghost = mp.has_longitudinal_field() ? add(fragments) : 0;

  double t = chooser.advance();
  for (std::vector<local_operator>::iterator opi = operators_p.begin();
       t < beta || opi != operators_p.end();) {

    // diagonal update & labeling
    if (opi == operators_p.end() || t < opi->time()) {
      // insert diagonal operator and graph if compatible
      local_graph g = chooser.diagonal();
      if ((is_bond(g) && is_compatible(g, spins_c[vsource(pos(g), vlattice())],
                                       spins_c[vtarget(pos(g), vlattice())])) ||
          (is_site(g) && is_compatible(g, spins_c[pos(g)]))) {
        operators.push_back(local_operator(g, t));
        t += chooser.advance();
      } else {
        t += chooser.advance();
        continue;
      }
    } else {
      if (opi->is_diagonal()) {
        // remove diagonal operator with probability one (= nothing to do)
        ++opi;
        continue;
      } else {
        // assign graph to offdiagonal operator
        opi->assign_graph(chooser.offdiagonal(opi->loc()));
        operators.push_back(*opi);
        ++opi;
      }
    }

    std::vector<local_operator>::reverse_iterator oi = operators.rbegin();
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
    for (boost::tie(rsi, rsi_end) = sites(rgraph());
         rsi != rsi_end; ++rsi) {
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
  for (std::vector<cluster_fragment>::iterator ci = fragments.begin();
       ci != fragments.end(); ++ci)
    if (ci->is_root()) {
      ci->id = clusters.size();
      clusters.push_back(cluster_info(random() < 0.5));
    }
  if (mp.has_longitudinal_field())
    clusters[cluster_id(fragments, ghost)].to_flip = false;

  // flip operators and spins & do improved measurements
  if (!mp.has_longitudinal_field()) {
    std::copy(spins.begin(), spins.end(), spins_c.begin());
    for (std::vector<local_operator>::iterator oi = operators.begin();
         oi != operators.end(); ++oi) {
      int id_l = root(fragments, oi->loop0).id;
      int id_u = root(fragments, oi->loop1).id;
      if (oi->is_bond()) {
        int s0 = vsource(oi->pos(), vlattice());
        int s1 = vtarget(oi->pos(), vlattice());
        // clusters[id_l].mag += 0 * oi->time;
        clusters[id_l].length += 2 * oi->time();
        if (oi->is_offdiagonal()) {
          spins[s0] ^= 1;
          spins[s1] ^= 1;
        }
        // clusters[id_u].mag -= 0 * oi->time;
        clusters[id_u].length -= 2 * oi->time();
      } else {
        int s = oi->pos();
        clusters[id_l].mag += (1 - 2 * spins[s]) * oi->time();
        clusters[id_l].length += oi->time();
        if (oi->is_offdiagonal()) spins[s] ^= 1;
        clusters[id_u].mag -= (1 - 2 * spins[s]) * oi->time();
        clusters[id_u].length -= oi->time();
      }
      if (clusters[id_l].to_flip ^ clusters[id_u].to_flip) oi->flip();
    }
    for (int s = 0; s < nvs; ++s) {
      int id = cluster_id(fragments, s);
      clusters[id].mag0 += (1 - 2 * spins[s]);
      clusters[id].size += 1;
      clusters[id].mag += (1 - 2 * spins[s]) * beta;
      clusters[id].length += beta;
      if (clusters[id].to_flip) spins[s] ^= 1;
    }
  } else {
    for (std::vector<local_operator>::iterator oi = operators.begin();
         oi != operators.end(); ++oi)
      if (clusters[cluster_id(fragments, oi->loop0)].to_flip ^
          clusters[cluster_id(fragments, oi->loop1)].to_flip) oi->flip();
    for (int s = 0; s < nvs; ++s)
      if (clusters[cluster_id(fragments, s)].to_flip) spins[s] ^= 1;
  }

  //
  // measurements
  //

  int nrs = num_sites(rgraph());

  // energy
  double energy = energy_offset;
  {
    bond_iterator rbi, rbi_end;
    for (boost::tie(rbi, rbi_end) = bonds(rgraph());
         rbi != rbi_end; ++rbi) {
      double jz = mp.bond(*rbi, rgraph()).jz;
      bond_iterator vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = virtual_bonds(vlattice(), rgraph(), *rbi);
           vbi != vbi_end; ++vbi) {
        if (spins_c[vsource(*vbi, vlattice())] == spins_c[vtarget(*vbi, vlattice())])
          energy += 0.25 * jz;
        else
          energy -= 0.25 * jz;
      }
      if (mp.has_d_term()) {
        site_iterator rsi, rsi_end;
        for (boost::tie(rsi, rsi_end) = sites(rgraph());
             rsi != rsi_end; ++rsi) {
          double d = mp.site(*rsi, rgraph()).d;
          bond_iterator vbi, vbi_end;
          for (boost::tie(vbi, vbi_end) =
                 virtual_bonds(vlattice(), rgraph(), *rsi);
               vbi != vbi_end; ++vbi) {
            if (spins_c[vsource(*vbi, vlattice())] == spins_c[vtarget(*vbi, vlattice())])
              energy += 0.5 * d;
            else
              energy -= 0.5 * d;
          }
        }
      }
    }
    int nop = 0;
    for (std::vector<local_operator>::iterator oi = operators.begin();
         oi != operators.end(); ++oi) if (oi->is_offdiagonal()) ++nop;
    energy -= (double)nop / beta;
  }
  measurements["Energy"] << energy;
  measurements["Energy Density"] << energy / nrs;

  {
    // accumurate loop length and magnetization
    double z2 = 0;
    double s2 = 0;
    double m2 = 0;
    double l2 = 0;
    std::vector<cluster_info>::iterator x;
    for (std::vector<cluster_info>::iterator pi = clusters.begin();
         pi != clusters.end(); ++pi) {
      z2 += looper::sqr(pi->mag0);
      s2 += looper::sqr(pi->size);
      m2 += looper::sqr(pi->mag);
      l2 += looper::sqr(pi->length);
    }

    measurements["Magnetization^2"] << z2 / (4 * nrs);
    measurements["Susceptibility"] << m2 / (4 * beta * nrs);
    if (is_bipartite()) {
      measurements["Staggered Magnetization^2"] << s2 / (4 * nrs);
      measurements["Staggered Susceptibility"]
        << l2 / (4 * beta * nrs);
    }
  }
}
