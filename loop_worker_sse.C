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

qmc_worker_sse::qmc_worker_sse(const alps::ProcessList& w,
                               const alps::Parameters& p, int n)
  : super_type(w, p, n),
    nrs(num_sites(rlat())), nvs(num_sites(vlat())),
    beta(1.0 / static_cast<double>(p["T"])),
    spins(nvs, 0 /* all up */), operators(0), nop(0),
    spins_curr(nvs), fragments(), current(nvs), clusters()
{}

void qmc_worker_sse::dostep()
{
  super_type::dostep();

  //
  // adjust length of operator string
  //

  if (nop > 0.8 * operators.size()) {
    std::vector<local_operator> operators_new;
    std::vector<local_operator>::iterator itr_new = operators_new.begin();
    for (std::vector<local_operator>::iterator itr = operators.begin();
         itr!= operators.end(); ++itr) {
      operators_new.push_back(*itr);
      operators_new.push_back(local_operator());
    }
    std::swap(operators, operators_new);
  }

  //
  // diagonal update and cluster construction
  //

  // initialize spin & operator information
  std::copy(spins.begin(), spins.end(), spins_curr.begin());

  // initialize cluster information (setup cluster fragments)
  fragments.resize(0); fragments.resize(nvs);
  for (int s = 0; s < nvs; ++s) current[s] = s;
  int ghost = has_longitudinal_field() ? add(fragments) : 0;

  for (std::vector<local_operator>::iterator oi = operators.begin();
       oi != operators.end(); ++oi) {

    // diagonal update & labeling
    if (oi->is_identity()) {
      // insert diagonal operator and graph if compatible
      local_graph g = choose_graph();
#warning "to be checked"
      if (((operators.size() - nop) * random() < beta / 2) &&
          ((is_bond(g) &&
            is_compatible(g, spins_curr[vsource(pos(g))],
                          spins_curr[vtarget(pos(g))])) ||
           (is_site(g) && is_compatible(g, spins_curr[pos(g)])))) {
        *oi = g;
        ++nop;
      } else
        continue;
    } else if (oi->is_diagonal()) {
      // remove diagonal operator with a certain probability
#warning "to be checked"
      if (beta / 2 * random() < operators.size() - nop + 1) {
        *oi = local_operator();
        --nop;
        continue;
      } else {
        // assign graph to diagonal operator
      }
    } else {
      // assign graph to offdiagonal operator
    }

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
  site_iterator rsi, rsi_end;
  for (boost::tie(rsi, rsi_end) = sites(rlat()); rsi != rsi_end; ++rsi) {
    site_iterator vsi, vsi_end;
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

}
