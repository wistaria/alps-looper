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

class bond_weight {

  // loop equations:
  //
  //   - offset + v1 + v2 = - C - Jz/4
  //   - offset + v3 + v4 = - C + Jz/4
  //              v1 + v3 =       |Jxy|/2

  // standard solution:
  //
  // i) Jz <= -|Jxy|  (ferro-Ising)
  //      v1 = |Jxy|/2
  //      v2 = -(|Jxy| + Jz)/2
  //      v3 = v4 = 0
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //      v1 = (|Jxy| - Jz)/4
  //      v3 = (|Jxy| + Jz)/4
  //      v2 = v4 = 0
  // iii) Jz >= |Jxy|  (antiferro-Ising)
  //      v3 = |Jxy|/2
  //      v4 = -(|Jxy| - Jz)/2
  //      v1 = v2 = 0

  // "ergodic" solutions (with additional parameter 0 < a < 1)
  //
  // i) Jz <= -|Jxy|  (ferro-Ising)
  //      same as the standard solution
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //   ii-1) |Jxy| - Jz >= 2 a |Jxy|
  //      same as the standard solution
  //   ii-2) |Jxy| - Jz <= 2 a |Jxy|
  //      v1 = a |Jxy| / 2
  //      v2 = 0
  //      v3 = (1-a) |Jxy| / 2
  //      v4 = -((1-2a) |Jxy| - Jz)/2
  // iii) Jz >= |Jxy|  (antiferro-Ising)
  //      same as ii-2)

public:
  bond_weight() : sign_(1), offset_(0) { v_[1] = v_[2] = v_[3] = v_[4] = 0; }
  bond_weight(const bond_parameter& p, double force_scatter = 0)
  { init(p, force_scatter); }

  void init(const bond_parameter& p, double force_scatter = 0)
  {
    sign_ = (p.jxy() <= 0 ? 1 : -1);
    double c = p.c();
    double jxy = std::abs(p.jxy());
    double jz = p.jz();
    double a = crop_01(force_scatter);
    if (alps::is_nonzero<1>(jxy + std::abs(jz))) {
      if (jxy - jz > 2 * a * jxy) {
        // standard solutions
        v_[1] = crop_0(std::min(jxy/2, (jxy - jz)/4));
        v_[2] = crop_0(-(jxy + jz)/2);
        v_[3] = crop_0(std::min(jxy/2, (jxy + jz)/4));
        v_[4] = crop_0(-(jxy - jz)/2);
      } else {
        // "ergodic" solutions
        v_[1] = a*jxy/2;
        v_[2] = 0;
        v_[3] = (1-a)*jxy/2;
        v_[4] = -((1-2*a)*jxy-jz)/2;
      }
    } else {
      v_[1] = v_[2] = v_[3] = v_[4] = 0;
    }
    offset_ = c + (v_[1] + v_[2] + v_[3] + v_[4])/2;
  }

  bool has_weight() const
  {
    return alps::is_nonzero<1>(v_[1] + v_[2] + v_[3] + v_[4]) &&
      (v_[1] + v_[2] + v_[3] + v_[4]) > 0;
  }
  double sign() const { return sign_; }
  double offset() const { return offset_; }
  double v(int g) const { return v_[g]; }

  static bond_parameter check(const bond_parameter& p, const bond_weight& w)
  {
    if (!alps::is_equal<1>(-w.offset()+w.v(1)+w.v(2), -p.c()-p.jz()/4) ||
        !alps::is_equal<1>(-w.offset()+w.v(3)+w.v(4), -p.c()+p.jz()/4) ||
        !alps::is_equal<1>(w.v(1)+w.v(3), -w.sign()*p.jxy()/2))
      boost::throw_exception(std::logic_error("bond_parameter::check() 1"));
    double c = w.offset() - (w.v(1) + w.v(2) + w.v(3) + w.v(4))/2;
    double jxy = -2 * (w.v(1) + w.v(3)) * w.sign();
    double jz = -2 * (w.v(1) + w.v(2) - w.v(3) - w.v(4));
    if (!alps::is_equal<1>(p.c(), c) || !alps::is_equal<1>(p.jxy(), jxy) ||
        !alps::is_equal<1>(p.jz(), jz)) {
      std::cerr << p.c() << ' ' << p.jxy() << ' ' << p.jz() << std::endl;
      std::cerr << c << ' ' << jxy << ' ' << jz << std::endl;
      boost::throw_exception(std::logic_error("bond_parameter::check() 2"));
    }
    return bond_parameter(c, jxy, jz);
  }

private:
  double sign_;
  double offset_;
  double v_[5]; // v_[0] is not used
};


class site_weight {

  // loop equations:
  //
  //   - offset + v1 + v2 = - C + Hz/2
  //   - offset + v1 + v3 = - C - Hz/2
  //              v1      =       |Hx|/2

  // standard solution:
  //
  // i) Hz >= 0
  //      v1 = |Hx|/2
  //      v2 = Hz/2
  //      v3 = 0
  // ii) Hz < 0
  //      v1 = |Hx|/2
  //      v2 = 0
  //      v3 = -Hz/2

public:
  site_weight() : sign_(1), offset_(0) { v_[1] = v_[2] = v_[3] = 0; }
  site_weight(const site_parameter& p) { init(p); }

  void init(const site_parameter& p)
  {
    sign_ = (p.hx() >= 0 ? 1 : -1);
    v_[1] = std::abs(p.hx()) / 2;
    v_[2] = crop_0( p.hz());
    v_[3] = crop_0(-p.hz());
    offset_ = p.c() + v_[1] + (v_[2] + v_[3])/2;
  }

  bool has_weight() const
  {
    return alps::is_nonzero<1>(v_[1] + v_[2] + v_[3]) &&
      (v_[1] + v_[2] + v_[3]) > 0;
  }
  double sign() const { return sign_; }
  double offset() const { return offset_; }
  double v(int g) const { return v_[g]; }

  static site_parameter check(const site_parameter& p, const site_weight& w)
  {
    if (!alps::is_equal<1>(-w.offset()+w.v(1)+w.v(2), -p.c()+p.hz()/2) ||
        !alps::is_equal<1>(-w.offset()+w.v(1)+w.v(3), -p.c()-p.hz()/2) ||
        !alps::is_equal<1>(w.v(1), w.sign()*p.hx()/2))
      boost::throw_exception(std::logic_error("site_parameter::check() 1"));
    double c = w.offset() - (w.v(1) + (w.v(2) + w.v(3))/2);
    double hx = 2 * w.v(1) * w.sign();
    double hz = w.v(2) - w.v(3);
    if (!alps::is_equal<1>(p.c(), c) || !alps::is_equal<1>(p.hx(), hx) ||
        !alps::is_equal<1>(p.hz(), hz)) {
      std::cerr << p.c() << ' ' << p.hx() << ' ' << p.hz() << std::endl;
      std::cerr << c << ' ' << hx << ' ' << hz << std::endl;
      boost::throw_exception(std::logic_error("site_parameter::check() 2"));
    }
    return site_parameter(0.5, c, hx, hz, 0);
  }

private:
  double sign_;
  double offset_;
  double v_[4]; // v_[0] is not used
};

} // namespace looper

#endif // LOOPER_WEIGHT_H
