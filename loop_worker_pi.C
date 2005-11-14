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

qmc_worker_pi::qmc_worker_pi(const alps::ProcessList& w,
			     const alps::Parameters& p, int n)
  : super_type(w, p, n),
    nrs(num_sites(rlat())), nvs(num_sites(vlat())), 
    beta(1.0 / static_cast<double>(p["T"])),
    spins(nvs, 0 /* all up */), operators(0), 
    spins_curr(nvs), operators_prev(), fragments(), current(nvs), clusters()
{}

void qmc_worker_pi::dostep()
{
  dostep();

  //
  // diagonal update and cluster construction
  //

  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_curr.begin());
  std::swap(operators, operators_prev); operators.resize(0);

  // initialize cluster information (setup cluster fragments)
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;
  int ghost = has_longitudinal_field() ? add(fragments) : 0;

  double t = advance();
  for (std::vector<local_operator>::iterator opi = operators_prev.begin();
       t < beta || opi != operators_prev.end();) {

    // diagonal update & labeling
    if (opi == operators_prev.end() || t < opi->time()) {
      // insert diagonal operator and graph if compatible
      looper::local_graph g = choose_graph();
      if ((is_bond(g) &&
           is_compatible(g, spins_curr[vsource(pos(g))],
                         spins_curr[vtarget(pos(g))])) ||
          (is_site(g) && is_compatible(g, spins_curr[pos(g)]))) {
        operators.push_back(local_operator(g, t));
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
        opi->assign_graph(choose_graph(opi->loc()));
        operators.push_back(*opi);
        ++opi;
      }
    }

    std::vector<local_operator>::reverse_iterator oi = operators.rbegin();
    if (oi->is_bond()) {
      int b = oi->pos();
      int s0 = vsource(b);
      int s1 = vtarget(b);
      boost::tie(current[s0], current[s1], oi->loop0, oi->loop1) =
        reconnect(fragments, oi->graph(), current[s0], current[s1]);
      if (oi->is_offdiagonal()) {
        spins_curr[s0] ^= 1;
        spins_curr[s1] ^= 1;
      }
    } else {
      int s = oi->pos();
      boost::tie(current[s], oi->loop0, oi->loop1) =
        reconnect(fragments, oi->graph(), current[s]);
      if (oi->is_locked()) unify(fragments, ghost, current[s]);
      if (oi->is_offdiagonal()) spins_curr[s] ^= 1;
    }
  }

  // connect bottom and top cluster fragments after random permutation
  std::vector<int> r, c0, c1;
  super_type::site_iterator rsi, rsi_end;
  for (boost::tie(rsi, rsi_end) = sites(rlat()); rsi != rsi_end; ++rsi) {
    super_type::site_iterator vsi, vsi_end;
    boost::tie(vsi, vsi_end) = virtual_sites(vlat(), rlat(), *rsi);
    int offset = *vsi;
    int s2 = *vsi_end - *vsi;
    r.resize(s2); c0.resize(s2); c1.resize(s2);
    for (int i = 0; i < s2; ++i) {
      r[i] = i;
      c0[i] = spins[offset+i];
      c1[i] = spins_curr[offset+i];
    }
    looper::restricted_random_shuffle(r.begin(), r.end(), c0.begin(),
				      c0.end(), c1.begin(), c1.end(),
				      *engine_ptr);
    for (int i = 0; i < s2; ++i)
      unify(fragments, offset+i, current[offset+r[i]]);
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
  if (has_longitudinal_field())
    clusters[cluster_id(fragments, ghost)].to_flip = false;

  // flip operators and spins & do improved measurements
  if (!has_longitudinal_field()) {
    std::copy(spins.begin(), spins.end(), spins_curr.begin());
    for (std::vector<local_operator>::iterator oi = operators.begin();
         oi != operators.end(); ++oi) {
      int id_l = root(fragments, oi->loop0).id;
      int id_u = root(fragments, oi->loop1).id;
      if (oi->is_site()) {
        int s = oi->pos();
        clusters[id_l].mag += (1 - 2 * spins[s]) * oi->time();
        clusters[id_l].length += oi->time();
        if (oi->is_offdiagonal()) spins[s] ^= 1;
        clusters[id_u].mag -= (1 - 2 * spins[s]) * oi->time();
        clusters[id_u].length -= oi->time();
      } else {
        int b = oi->pos();
        int s0 = vsource(b);
        int s1 = vtarget(b);
        // clusters[id_l].mag += 0 * oi->time;
        clusters[id_l].length += 2 * oi->time();
        if (oi->is_offdiagonal()) {
          spins[s0] ^= 1;
          spins[s1] ^= 1;
        }
        // clusters[id_u].mag -= 0 * oi->time;
        clusters[id_u].length -= 2 * oi->time();
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
          clusters[cluster_id(fragments, oi->loop1)].to_flip)
        oi->flip();
    for (int s = 0; s < nvs; ++s)
      if (clusters[cluster_id(fragments, s)].to_flip) spins[s] ^= 1;
  }

  //
  // measurements
  //

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

    measurements["Energy"]
      << - (double)operators.size() / beta / nrs;
    measurements["Uniform Magnetization^2"] << z2 / (4 * nrs);
    measurements["Staggered Magnetization^2"] << s2 / (4 * nrs);
    measurements["Uniform Susceptibility"] << m2 / (4 * beta * nrs);
    measurements["Staggered Susceptibility"]
      << l2 / (4 * beta * nrs);
  }
}
