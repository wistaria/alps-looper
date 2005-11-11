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
  { return alps::is_nonzero<1>(weight()) && weight() > 0; }

  void check(const site_parameter& p) const;
};

struct bond_weight {

  double sign;
  double offset;
  double v[0];

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
        v[1] = crop_0(-(jxy - jz)/2);
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
  { return alps::is_nonzero<1>(weight()) && weight() > 0; }

  void check(const bond_parameter& p) const;
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
