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
// #include <alps/fixed_capacity_vector.h>

qmc_worker_sse::qmc_worker_sse(const alps::ProcessList& w,
                               const alps::Parameters& p, int n)
  : super_type(w, p, n) /* , chooser(*engine_ptr) */
{
//   beta = 1.0 / static_cast<double>(p["T"]);
//   if (beta < 0)
//     boost::throw_exception(std::invalid_argument("negative beta"));

//   //
//   // setup model
//   //

//   looper::model_parameter mp(p, *this);
//   is_signed = mp.is_signed();
//   is_classically_frustrated = mp.is_classically_frustrated();
//   has_hz = mp.has_longitudinal_field();
//   if (is_signed)
//     std::cerr << "WARNING: model has negative signs\n";
//   if (is_classically_frustrated)
//     std::cerr << "WARNING: model is classically frustrated\n";

//   energy_offset = 0;
//   {
//     alps::graph_traits<graph_type>::site_iterator si, si_end;
//     for (boost::tie(si, si_end) = sites(real_graph()); si != si_end; ++si)
//       energy_offset += mp.site(*si, real_graph()).c;
//     alps::graph_traits<graph_type>::bond_iterator bi, bi_end;
//     for (boost::tie(bi, bi_end) = bonds(real_graph()); bi != bi_end; ++bi)
//       energy_offset += mp.bond(*bi, real_graph()).c;
//   }

//   //
//   // setup virtual lattice
//   //

//   is_bipartite = alps::set_parity(real_graph());
//   vlat.generate(real_graph(), mp, mp.has_d_term());

//   //
//   // setup graph chooser
//   //

//   chooser.init(looper::weight_table(mp, real_graph(), vlat));

//   //
//   // initialize configuration
//   //

//   int nvs = num_sites(vlat);
//   spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);
//   operators.resize(0);
//   nop = 0;
//   spins_c.resize(nvs);
//   current.resize(nvs);

//   //
//   // init measurements
//   //

//   using alps::RealObservable;
//   using alps::make_observable;

//   if (is_signed) {
//     measurements << RealObservable("Sign");
//   }

//   measurements
//     << make_observable(
//          RealObservable("Energy"), is_signed)
//     << make_observable(
//          RealObservable("Energy Density"), is_signed)
//     << make_observable(
//          RealObservable("Diagonal Energy Density"), is_signed)
//     << make_observable(
//          RealObservable("Energy Density^2"), is_signed)
//     << make_observable(
//          RealObservable("beta * Energy / sqrt(N)"), is_signed)
//     << make_observable(
//          RealObservable("beta * Energy^2"), is_signed);

//   measurements
//     << make_observable(
//          RealObservable("Magnetization"), is_signed)
//     << make_observable(
//          RealObservable("Magnetization^2"), is_signed)
//     << make_observable(
//          RealObservable("Susceptibility"), is_signed);

//   if (is_bipartite) {
//     measurements
//       << make_observable(
//            RealObservable("Staggered Magnetization"), is_signed)
//       << make_observable(
//            RealObservable("Staggered Magnetization^2"), is_signed)
//       << make_observable(
//            RealObservable("Staggered Susceptibility"), is_signed);
//   }

//   if (!is_classically_frustrated) {
//     measurements
//       << RealObservable("Generalized Magnetization^2")
//       << RealObservable("Generalized Susceptibility");
//   }
//   if (!is_classically_frustrated && is_bipartite) {
//     measurements
//       << RealObservable("Staggered Generalized Magnetization^2")
//       << RealObservable("Staggered Generalized Susceptibility");
//   }
}

void qmc_worker_sse::dostep()
{
//   super_type::dostep();

//   //
//   // adjust length of operator string
//   //

//   if (nop > 0.8 * operators.size()) {
//     std::vector<local_operator> operators_new;
//     std::vector<local_operator>::iterator itr_new = operators_new.begin();
//     for (std::vector<local_operator>::iterator itr = operators.begin();
//          itr!= operators.end(); ++itr) {
//       operators_new.push_back(*itr);
//       operators_new.push_back(local_operator());
//     }
//     std::swap(operators, operators_new);
//   }

//   //
//   // diagonal update and cluster construction
//   //

//   // initialize spin & operator information
//   std::copy(spins.begin(), spins.end(), spins_c.begin());

//   // initialize cluster information (setup cluster fragments)
//   int nvs = num_sites(vlat);
//   fragments.resize(0); fragments.resize(nvs);
//   for (int s = 0; s < nvs; ++s) current[s] = s;
//   int ghost = has_hz ? add(fragments) : 0;

//   for (std::vector<local_operator>::iterator oi = operators.begin();
//        oi != operators.end(); ++oi) {

//     // diagonal update & labeling
//     if (oi->is_identity()) {
//       // insert diagonal operator and graph if compatible
//       local_graph g = chooser.diagonal();
// #warning "to be checked"
//       if (((operators.size() - nop) * random() < beta / 2) &&
//           ((is_bond(g) &&
//             is_compatible(g, spins_c[vsource(pos(g), vlat)],
//                           spins_c[vtarget(pos(g), vlat)])) ||
//            (is_site(g) && is_compatible(g, spins_c[pos(g)])))) {
//         *oi = g;
//         ++nop;
//       } else
//         continue;
//     } else if (oi->is_diagonal()) {
//       // remove diagonal operator with a certain probability
// #warning "to be checked"
//       if (beta / 2 * random() < operators.size() - nop + 1) {
//         *oi = local_operator();
//         --nop;
//         continue;
//       } else {
//         // assign graph to diagonal operator
//       }
//     } else {
//       // assign graph to offdiagonal operator
//     }

//     if (oi->is_bond()) {
//       int s0 = vsource(oi->pos(), vlat);
//       int s1 = vtarget(oi->pos(), vlat);
//       boost::tie(current[s0], current[s1], oi->loop0, oi->loop1) =
//         reconnect(fragments, oi->graph(), current[s0], current[s1]);
//       if (oi->is_offdiagonal()) {
//         spins_c[s0] ^= 1;
//         spins_c[s1] ^= 1;
//       }
//     } else {
//       int s = oi->pos();
//       boost::tie(current[s], oi->loop0, oi->loop1) =
//         reconnect(fragments, oi->graph(), current[s]);
//       if (oi->is_locked()) unify(fragments, ghost, current[s]);
//       if (oi->is_offdiagonal()) spins_c[s] ^= 1;
//     }
//   }

//   // connect bottom and top cluster fragments after random permutation
//   {
//     alps::fixed_capacity_vector<int, 20> r;
//     site_iterator rsi, rsi_end;
//     for (boost::tie(rsi, rsi_end) = sites(real_graph());
//          rsi != rsi_end; ++rsi) {
//       site_iterator vsi, vsi_end;
//       boost::tie(vsi, vsi_end) = virtual_sites(vlat, real_graph(), *rsi);
//       int offset = *vsi;
//       int s2 = *vsi_end - *vsi;
//       if (s2 == 1) {
//         unify(fragments, offset, current[offset]);
//       } else if (s2 > 1) {
//         r.resize(s2);
//         for (int i = 0; i < s2; ++i) r[i] = i;
//         looper::restricted_random_shuffle(r.begin(), r.end(),
//           spins.begin() + offset, spins_c.begin() + offset, random);
//         for (int i = 0; i < s2; ++i)
//           unify(fragments, offset+i, current[offset+r[i]]);
//       }
//     }
//   }

//   //
//   // cluster flip
//   //

//   // assign cluster id & determine if clusters are to be flipped
//   clusters.resize(0);
//   for (std::vector<cluster_fragment>::iterator ci = fragments.begin();
//        ci != fragments.end(); ++ci)
//     if (ci->is_root()) {
//       ci->id = clusters.size();
//       clusters.push_back(cluster_info(random() < 0.5));
//     }
//   if (has_hz)
//     clusters[cluster_id(fragments, ghost)].to_flip = false;

}
