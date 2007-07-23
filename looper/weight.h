/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#include "crop.h"
#include "model_parameter.h"

#include <alps/math.hpp>
#include <algorithm> // for std::min std::max
#include <cmath>     // for std::abs

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

struct site_weight {

  int sign;
  double offset;
  double v[1];

  // loop equations:
  //
  //   - offset + v0 = 0
  //   - offset + v0 = 0
  //              v0 = |Hx|/2

  // standard solution:
  //
  //   v0 = |Hx|/2

  site_weight() : sign(1), offset(0) { v[0] = 0; }
  site_weight(const site_parameter& p, double /* force_scatter */ = 0)
  { init(p); }

  void init(const site_parameter& p, double /* force_scatter */ = 0)
  {
    sign = (p.hx >= 0 ? 1 : -1);
    v[0] = std::abs(p.hx) / 2;
    offset = v[0];
  }

  double weight() const { return v[0]; }
  bool has_weight() const { return alps::is_positive<1>(weight()); }

  void check(const site_parameter& p) const;
};

struct bond_weight {

  int sign;
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
  //   - offset + v2 + v3 = - Jz/4
  //   - offset + v0 + v1 = + Jz/4
  //              v0 + v2 =   |Jxy|/2

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
  bond_weight(const site_parameter& p, double force_scatter = 0) { init(p, force_scatter); }
  bond_weight(const bond_parameter_xxz& p, double force_scatter = 0) { init(p, force_scatter); }
  bond_weight(const bond_parameter_xyz& p, double force_scatter = 0) { init(p, force_scatter); }

  void init(const bond_parameter_xxz& p, double force_scatter = 0) {
    sign = (p.jxy <= 0 ? 1 : -1);
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
    offset = weight()/2;
  }
  void init(const bond_parameter_xyz& p, double force_scatter = 0) {
    if (!alps::is_equal<>(p.jx, p.jy))
      boost::throw_exception(std::runtime_error("not an XXZ model"));
    init(bond_parameter_xxz(p.c, p.jx, p.jz), force_scatter);
  }
  void init(const site_parameter& p, double force_scatter = 0) {
    bond_parameter_xyz bp(0, 0, 0, p.d);
    init(bp, force_scatter);
  }

  double weight() const { return v[0] + v[1] + v[2] + v[3]; }
  bool has_weight() const
  { return alps::is_positive<1>(weight()); }

  void check(const bond_parameter_xxz& p) const;
};

class weight_table
{
public:
  typedef std::pair<int, site_weight> site_weight_t;
  typedef std::pair<int, bond_weight> bond_weight_t;
  typedef std::vector<std::pair<int, site_weight> >::const_iterator
    site_weight_iterator;
  typedef std::vector<std::pair<int, bond_weight> >::const_iterator
    bond_weight_iterator;

  template<class M, class LAT>
  weight_table(M const& m, const LAT& lat, double force_scatter = 0)
  { init(m, lat, force_scatter); }

  template<class M, class LAT>
  void init(M const& m, const LAT& lat, double force_scatter = 0)
  {
    weight_ = 0;
    site_weights_.clear();
    bond_weights_.clear();
    typename real_site_iterator<LAT>::type rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(lat.rg());
         rsi != rsi_end; ++rsi) {
      site_weight sw(m.site(*rsi, lat.rg()), force_scatter);
      typename virtual_site_iterator<LAT>::type vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) = sites(lat, *rsi);
           vsi != vsi_end; ++vsi) {
        site_weights_.push_back(
          std::make_pair(get(site_index_t(), lat.vg(), *vsi), sw));
        weight_ += sw.weight();
      }
    }
    typename real_bond_iterator<LAT>::type rbi, rbi_end;
    for (boost::tie(rbi, rbi_end) = bonds(lat.rg());
         rbi != rbi_end; ++rbi) {
      bond_weight bw(m.bond(*rbi, lat.rg()), force_scatter);
      typename virtual_bond_iterator<LAT>::type vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = bonds(lat, *rbi);
           vbi != vbi_end; ++vbi) {
        bond_weights_.push_back(
          std::make_pair(get(bond_index_t(), lat.vg(), *vbi), bw));
        weight_ += bw.weight();
      }
    }
    if (m.has_d_term()) {
      for (boost::tie(rsi, rsi_end) = sites(lat.rg());
           rsi != rsi_end; ++rsi) {
        bond_weight bw(m.site(*rsi, lat.rg()), force_scatter);
        typename virtual_bond_iterator<LAT>::type vbi, vbi_end;
        for (boost::tie(vbi, vbi_end) = bonds(lat, *rsi);
             vbi != vbi_end; ++vbi) {
          bond_weights_.push_back(
             std::make_pair(get(bond_index_t(), lat.vg(), *vbi), bw));
          weight_ += bw.weight();
        }
      }
    }

    energy_offset_ = 0;
    for (site_weight_iterator itr = site_weights_.begin();
         itr != site_weights_.end(); ++itr)
      energy_offset_ += itr->second.offset;
    for (bond_weight_iterator itr = bond_weights_.begin();
         itr != bond_weights_.end(); ++itr)
      energy_offset_ += itr->second.offset;
  }

  std::pair<site_weight_iterator, site_weight_iterator>
  site_weights() const
  { return std::make_pair(site_weights_.begin(), site_weights_.end()); }
  std::pair<bond_weight_iterator, bond_weight_iterator>
  bond_weights() const
  { return std::make_pair(bond_weights_.begin(), bond_weights_.end()); }

  double weight() const { return weight_; }
  double rho() const { return weight(); }
  double energy_offset() const { return energy_offset_; }

private:
  double weight_;
  double energy_offset_;
  std::vector<std::pair<int, site_weight> > site_weights_;
  std::vector<std::pair<int, bond_weight> > bond_weights_;
};

//
// Implementations
//

inline void site_weight::check(const site_parameter& p) const
{
  if (!alps::is_zero<1>(-offset+v[0]) ||
      !alps::is_zero<1>(-offset+v[0]) ||
      !alps::is_equal<1>(v[0], sign * p.hx/2))
    boost::throw_exception(std::logic_error("site_parameter::check 1"));
  site_parameter pp(p.s, p.c, 2 * v[0] * sign, 0, 0);
  if (pp != p) {
    std::cerr << p.c << ' ' << p.hx << ' ' << pp.hz << std::endl;
    std::cerr << pp.c << ' ' << pp.hx << ' ' << pp.hz << std::endl;
    boost::throw_exception(std::logic_error("site_parameter::check 2"));
  }
}

inline void bond_weight::check(const bond_parameter_xxz& p) const
{
  if (!alps::is_equal<1>(-offset+v[0]+v[1],  p.jz/4) ||
      !alps::is_equal<1>(-offset+v[2]+v[3], -p.jz/4) ||
      !alps::is_equal<1>(v[0]+v[2], -sign * p.jxy/2))
    boost::throw_exception(std::logic_error("bond_parameter_xxz::check 1"));
  bond_parameter_xxz pp(p.c, -2 * (v[0] + v[2]) * sign, 2 * (v[0] + v[1] - v[2] - v[3]));
  if (pp != p) {
    std::cerr << p.c << ' ' << p.jxy << ' ' << p.jz << std::endl;
    std::cerr << pp.c << ' ' << pp.jxy << ' ' << pp.jz << std::endl;
    boost::throw_exception(std::logic_error("bond_parameter_xxz::check 2"));
  }
}

} // namespace looper

#endif // LOOPER_WEIGHT_H
