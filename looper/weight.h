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

#ifndef LOOPER_WEIGHT_H
#define LOOPER_WEIGHT_H

#include <looper/random_choice.h>
#include <looper/util.h>
#include <alps/math.hpp>
#include <algorithm> // for std::min std::max
#include <cmath>     // for std::abs

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

struct site_weight {

  double sign;
  double offset;
  double v[3];

  // loop equations:
  //
  //   - offset + v0 + v1 = - C + Hz/2
  //   - offset + v0 + g2 = - C - Hz/2
  //              v0      =       |Hx|/2

  // standard solution:
  //
  // i) Hz >= 0
  //      v0 = |Hx|/2
  //      v1 = Hz/2
  //      v2 = 0
  // ii) Hz < 0
  //      v0 = |Hx|/2
  //      v1 = 0
  //      v2 = -Hz/2

  site_weight() : sign(1), offset(0) { v[0] = v[1] = v[2] = 0; }
  site_weight(const site_parameter& p) { init(p); }

  void init(const site_parameter& p)
  {
    sign = (p.hx >= 0 ? 1 : -1);
    v[0] = std::abs(p.hx) / 2;
    v[1] = crop_0( p.hz);
    v[2] = crop_0(-p.hz);
    offset = p.c + v[0] + (v[1] + v[2])/2;
  }

  double weight() const { return v[0] + v[1] + v[2]; }
  bool has_weight() const
  { return alps::is_positive<1>(weight()); }

  void check(const site_parameter& p) const;
};

struct bond_weight {

  double sign;
  double offset;
  double v[4];

  // correspondence with notation in textbook
  //        textbook
  // v0     g3
  // v1     g4
  // v2     g1
  // v3     g2

  // loop equations:
  //
  //   - offset + v2 + v3 = - C - Jz/4
  //   - offset + v0 + v1 = - C + Jz/4
  //              v0 + v2 =       |Jxy|/2

  // standard solution:
  //
  // i) Jz <= -|Jxy|  (ferro-Ising)
  //      v2 = |Jxy|/2
  //      v3 = -(|Jxy| + Jz)/2
  //      v0 = v1 = 0
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //      v0 = (|Jxy| + Jz)/4
  //      v2 = (|Jxy| - Jz)/4
  //      v1 = v3 = 0
  // iii) Jz >= |Jxy|  (antiferro-Ising)
  //      v0 = |Jxy|/2
  //      v1 = -(|Jxy| - Jz)/2
  //      v2 = v3 = 0

  // "ergodic" solutions (with additional parameter 0 < a < 1)
  //
  // i) Jz <= -|Jxy|  (ferro-Ising)
  //      same as the standard solution
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //   ii-1) |Jxy| - Jz >= 2 a |Jxy|
  //      same as the standard solution
  //   ii-2) |Jxy| - Jz <= 2 a |Jxy|
  //      v0 = (1-a) |Jxy| / 2
  //      v1 = -((1-2a) |Jxy| - Jz)/2
  //      v2 = a |Jxy| / 2
  //      v3 = 0
  // iii) Jz >= |Jxy|  (antiferro-Ising)
  //      same as ii-2)

  bond_weight() : sign(1), offset(0) { v[0] = v[1] = v[2] = v[3] = 0; }
  bond_weight(const site_parameter& p, double force_scatter = 0)
  { init(p, force_scatter); }
  bond_weight(const bond_parameter& p, double force_scatter = 0)
  { init(p, force_scatter); }

  void init(const bond_parameter& p, double force_scatter = 0)
  {
    sign = (p.jxy <= 0 ? 1 : -1);
    double c = p.c;
    double jxy = std::abs(p.jxy);
    double jz = p.jz;
    double a = crop_01(force_scatter);
    if (alps::is_nonzero<1>(jxy + std::abs(jz))) {
      if (jxy - jz > 2 * a * jxy) {
        // standard solutions
        v[0] = crop_0(std::min(jxy/2, (jxy + jz)/4));
        v[1] = crop_0(-(jxy - jz)/2.0);
        v[2] = crop_0(std::min(jxy/2, (jxy - jz)/4));
        v[3] = crop_0(-(jxy + jz)/2);
      } else {
        // "ergodic" solutions
        v[0] = (1-a)*jxy/2;
        v[1] = -((1-2*a)*jxy-jz)/2;
        v[2] = a*jxy/2;
        v[3] = 0;
      }
    } else {
      v[0] = v[1] = v[2] = v[3] = 0;
    }
    offset = c + weight()/2;
  }
  void init(const site_parameter& p, double force_scatter = 0)
  {
    bond_parameter bp(0, 0, p.d);
    init(bp, force_scatter);
  }

  double weight() const { return v[0] + v[1] + v[2] + v[3]; }
  bool has_weight() const
  { return alps::is_positive<1>(weight()); }

  void check(const bond_parameter& p) const;
};

class weight_table
{
public:
  template<class RL, class VL>
  weight_table(const model_parameter& mp, const RL& rl, const VL& vl)
  { init(mp, rl, vl); }

  template<class RL, class VL>
  void init(const model_parameter& mp, const RL& rl, const VL& vl)
  {
    weight_ = 0;
    site_weights_.clear();
    bond_weights_.clear();
    typename graph_traits<RL>::site_iterator rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(rl); rsi != rsi_end; ++rsi) {
      site_weight sw(mp.site(*rsi, rl));
      typename graph_traits<VL>::site_iterator vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) = virtual_sites(vl, rl, *rsi);
           vsi != vsi_end; ++vsi) {
        site_weights_.push_back(
          std::make_pair(get(site_index_t(), vl.graph(), *vsi), sw));
        weight_ += sw.weight();
      }
    }
    typename graph_traits<RL>::bond_iterator rbi, rbi_end;
    for (boost::tie(rbi, rbi_end) = bonds(rl); rbi != rbi_end; ++rbi) {
      bond_weight bw(mp.bond(*rbi, rl));
      typename graph_traits<VL>::bond_iterator vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = virtual_bonds(vl, rl, *rbi);
           vbi != vbi_end; ++vbi) {
        bond_weights_.push_back(
          std::make_pair(get(bond_index_t(), vl.graph(), *vbi), bw));
        weight_ += bw.weight();
      }
    }
    if (mp.has_d_term()) {
      for (boost::tie(rsi, rsi_end) = sites(rl); rsi != rsi_end; ++rsi) {
        bond_weight bw(mp.site(*rsi, rl));
        typename graph_traits<VL>::bond_iterator vbi, vbi_end;
        for (boost::tie(vbi, vbi_end) = virtual_bonds(vl, rl, *rsi);
             vbi != vbi_end; ++vbi) {
          bond_weights_.push_back(
             std::make_pair(get(bond_index_t(), vl.graph(), *vbi), bw));
          weight_ += bw.weight();
        }
      }
    }
  }

  double weight() const { return weight_; }
  double rho() const { return weight(); }

  template<class LOCAL_GRAPH, class RND_CHOICE>
  void setup_graph_chooser(std::vector<LOCAL_GRAPH>& diagonal_graphs,
                           RND_CHOICE& rng,
                           std::vector<double>& offdiagonal_weights) const
  {
    diagonal_graphs.clear();
    offdiagonal_weights.clear();
    std::vector<double> w;
    for (std::vector<std::pair<int, site_weight> >::const_iterator
           itr = site_weights_.begin(); itr != site_weights_.end(); ++itr) {
      site_weight sw = itr->second;
      for (int g = 0; g <= 2; ++g) {
        if (alps::is_nonzero<1>(sw.v[g])) {
          diagonal_graphs.push_back(LOCAL_GRAPH::site_graph(g, itr->first));
          w.push_back(sw.v[g]);
        }
      }
    }
    for (std::vector<std::pair<int, bond_weight> >::const_iterator
           itr = bond_weights_.begin(); itr != bond_weights_.end(); ++itr) {
      bond_weight bw = itr->second;
      for (int g = 0; g <= 3; ++g) {
        if (alps::is_nonzero<1>(bw.v[g])) {
          diagonal_graphs.push_back(LOCAL_GRAPH::bond_graph(g, itr->first));
          w.push_back(bw.v[g]);
        }
      }
      offdiagonal_weights.push_back(alps::is_nonzero<1>(bw.v[0] + bw.v[2]) ?
                                   bw.v[0] / (bw.v[0] + bw.v[2]) : 1);
    }
    rng.distribution().init(w);
  }

private:
  double weight_;
  std::vector<std::pair<int, site_weight> > site_weights_;
  std::vector<std::pair<int, bond_weight> > bond_weights_;
};

//
// Implementations
//

inline void site_weight::check(const site_parameter& p) const
{
  if (!alps::is_equal<1>(-offset+v[0]+v[1], -p.c+p.hz/2) ||
      !alps::is_equal<1>(-offset+v[0]+v[2], -p.c-p.hz/2) ||
      !alps::is_equal<1>(v[0], sign*p.hx/2))
    boost::throw_exception(std::logic_error("site_parameter::check 1"));
  site_parameter pp(p.s,
                    offset - (v[0] + (v[1] + v[2])/2),
                    2 * v[0] * sign,
                    v[1] - v[2]);
  if (pp != p) {
    std::cerr << p.c << ' ' << p.hx << ' ' << p.hz << std::endl;
    std::cerr << pp.c << ' ' << pp.hx << ' ' << pp.hz << std::endl;
    boost::throw_exception(std::logic_error("site_parameter::check 2"));
  }
}

inline void bond_weight::check(const bond_parameter& p) const
{
  if (!alps::is_equal<1>(-offset+v[0]+v[1], -p.c+p.jz/4) ||
      !alps::is_equal<1>(-offset+v[2]+v[3], -p.c-p.jz/4) ||
      !alps::is_equal<1>(v[0]+v[2], -sign*p.jxy/2))
    boost::throw_exception(std::logic_error("bond_parameter::check 1"));
  bond_parameter pp(offset - weight()/2,
                    -2 * (v[0] + v[2]) * sign,
                    2 * (v[0] + v[1] - v[2] - v[3]));
  if (pp != p) {
    std::cerr << p.c << ' ' << p.jxy << ' ' << p.jz << std::endl;
    std::cerr << pp.c << ' ' << pp.jxy << ' ' << pp.jz << std::endl;
    boost::throw_exception(std::logic_error("bond_parameter::check 2"));
  }
}

} // namespace looper

#endif // LOOPER_WEIGHT_H
