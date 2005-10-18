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
  //              v1 + v3 = |Jxy|/2

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
  bond_weight() : offset_(0), sign_(1)
  { v_[1] = v_[2] = v_[3] = v_[4] = 0; }
  bond_weight(const bond_parameter& p, double force_scatter = 0)
  { init(p, force_scatter); }

  void init(const bond_parameter& p, double force_scatter = 0)
  {
    using alps::is_nonzero;
    sign_ = (p.jxy() >= 0 ? 1 : -1);
    double c = p.c();
    double jxy = std::abs(p.jxy());
    double jz = p.jz();
    double a = crop_01(force_scatter);
    if (is_nonzero<1>(jxy + std::abs(jz))) {
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
        v_[4] = -((1-2*z)*jxy-jz)/2;
      }
    } else {
      v_[1] = v_[2] = v_[3] = v_[4] = 0;
    }
    offset_ = c + (v_[1] + v_[2] + v_[3] + v_[4])/2;
  }

  bool weight() const
  { return is_nonzero<1>(v_[1] + v_[2] + v_[3] + v_[4]); }
  double v(int g) const { return v_[g]; }
  double offset() const { return offset_; }
  double sign() const { return sign_; }

  static bond_parameter check(const bond_parameter& p, const bond_weight& w)
  {
    using alps::is_equal;
    if (!is_equal<1>(-w.offset() + w.v(1) + w.v(2) = -p.c() - p.jz()/4) || 
	!is_equal<1>(-w.offset() + w.v(3) + w.v(4) = -p.c() + p.jz()/4) || 
	!is_equal<1>(w.v(1) + w.v(3) = - w.sign() * p.jxy()/4))
      boost::throw_exception(std::logic_error());
    double c = w.offset() - (w.v(1) + w.v(2) + w.v(3) + w.v(4));
    double jxy = -2 * (w.v(1) + w.v(3)) * w.sign();
    double jz = -2 * (w.v(1) + w.v(2) - w.v(3) - w.v(4));
    if (!is_equal<1>(p.c(), c) || !is_equal<1>(p.jz(), jz) ||
	!is_equal<1>(p.jxy(), jxy))
      boost::throw_exception(std::logic_error());
    return bond_parameter(c, jxy, jz);
  }

private:
  double offset_;
  double sign_;
  double v_[5]; // v_[0] is not used
};


class site_weight {

  // loop equation and its solution:
  //
  //   w1 = Hx/2
  //
  //   density = w1

public:
  site_weight() : v_(0), sign_(1) {}
  site_weight(const site_parameter& p) { init(p); }

  void init(const site_parameter& p)
  {
    v_ = std::abs(p.hx()) / 2;
    sign_ = (p.hx() >= 0 ? 1 : -1);
  }

  double weight() const { return v_; }
  double offset() const { return 0; }
  double sign() const { return sign_; }

  static site_parameter check(const site_weight& w)
  {
    double hx = 2 * w.weight() * w.sign();
    return site_parameter(0, hx);
  }

private:
  double v_;
  double sign_;
};


template<class WEIGHT>
class uniform_bond_chooser
{
public:
  typedef WEIGHT weight_type;

  uniform_bond_chooser() : n_(), weight_(), gw_() {}
  template<class GRAPH, class MODEL>
  uniform_bond_chooser(const GRAPH& g, const MODEL& m, double fs = 0)
    : n_(), weight_(), gw_()
  { init(g, m, fs); }

  template<class GRAPH, class MODEL>
  void init(const GRAPH& g, const MODEL& m, double fs = 0)
  {
    assert(m.num_bond_types() == 1);
    n_ = double(boost::num_edges(g));
    weight_ = weight_type(m.uniform_bond(), fs);
    gw_ = n_ * weight_.weight();
  }

  template<class RNG>
  int choose(RNG& rng) const { return n_ * rng(); }
  template<class RNG>
  int operator()(RNG& rng) const { return choose(rng); }

  weight_type& weight(int) { return weight_; }
  const weight_type& weight(int) const { return weight_; }
  double global_weight() const { return gw_; }

private:
  double n_;
  weight_type weight_;
  double gw_;
};

template<class WEIGHT>
class bond_chooser
{
public:
  typedef WEIGHT weight_type;

  bond_chooser() : weight_(), rc_(), gw_(0) {}
  template<class G, class M>
  bond_chooser(const G& rg, const G& vg, const virtual_mapping<G>& vm,
               const M& m, double fs = 0)
    : weight_(), rc_(), gw_(0)
  { init (rg, vg, vm, m, fs); }

  template<class G, class M>
  void init(const G& rg, const G& vg, const virtual_mapping<G>& vm,
            const M& m, double fs = 0)
  {
    weight_.clear();
    gw_ = 0.0;
    if (boost::num_edges(vg) > 0) {
      typename boost::graph_traits<G>::edge_iterator rei, rei_end;
      for (boost::tie(rei, rei_end) = boost::edges(rg);
           rei != rei_end; ++rei) {
        typename boost::graph_traits<G>::edge_iterator vei, vei_end;
        for (boost::tie(vei, vei_end) = vm.virtual_edges(rg, *rei);
             vei != vei_end; ++vei) {
          weight_.push_back(weight_type(m.bond(*rei, rg), fs));
        }
      }

      std::vector<double> w(0);
      typename std::vector<weight_type>::iterator itr_end = weight_.end();
      for (typename std::vector<weight_type>::iterator itr = weight_.begin();
           itr != itr_end; ++itr) {
        w.push_back(itr->weight());
        gw_ += itr->weight();
      }
      rc_.init(w);
    }
  }

  template<class RNG>
  int choose(RNG& rng) const { return rc_(rng); }
  template<class RNG>
  int operator()(RNG& rng) const { return choose(rng); }

  weight_type& weight(int i) { return weight_[i]; }
  const weight_type& weight(int i) const { return weight_[i]; }
  double global_weight() const { return gw_; }

private:
  std::vector<weight_type> weight_;
  random_choice<> rc_;
  double gw_;
};

} // namespace looper

#endif // LOOPER_WEIGHT_H
