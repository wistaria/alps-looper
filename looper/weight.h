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
  //   w0 + w1 + w2 = Jz/4
  //   w0 + w3 + w4 = -Jz/4
  //   w1 + w3 = |Jxy|/2
  //
  //   density = max((w1 + w2), (w3 + w4))
  //   Pa_para = (w1 + w2) / density
  //   Pa_anti = (w3 + w4) / density
  //   Pf_para = w2 / (w1 + w2)
  //   Pf_anti = w4 / (w3 + w4)
  //   Pr = w3 / (w1 + w3)

  // standard solution:
  //
  // i) Jz >= |Jxy|  (ferro-Ising)
  //      w1 = |Jxy|/2
  //      w2 = (Jz - |Jxy|)/2
  //      w3 = w4 = 0
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //      w1 = (Jz + |Jxy|)/4
  //      w3 = (-Jz + |Jxy|)/4
  //      w2 = w4 = 0
  // iii) Jz <= -|Jxy|  (antiferro-Ising)
  //      w3 = |Jxy|/2
  //      w4 = (-Jz - |Jxy|)/2
  //      w1 = w2 = 0

  // "ergodic" solutions (with additional parameter a)
  //
  // i) Jz >= |Jxy|  (ferro-Ising)
  //      same as the standard solution
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //   ii-1) Jz + |Jxy| >= 2 a |Jxy|
  //      same as the standard solution
  //   ii-2) Jz + |Jxy| <= 2 a |Jxy|
  //      w1 = a |Jxy| / 2
  //      w2 = 0
  //      w3 = (1-a) |Jxy| / 2
  //      w4 = (-Jz + (2a-1) |Jxy|)/2
  // iii) Jz <= -|Jxy|  (antiferro-Ising)
  //      same as ii-2)

public:
  bond_weight() : weight_(0), offset_(0), sign_(1) 
  { v[1] = v[2] = v[3] = v[4] = 0; }
  bond_weight(const bond_parameter& p, double force_scatter = 0) 
  { init(p, force_scatter); }

  void init(const bond_parameter& p, double force_scatter = 0)
  {
    double Jxy = std::abs(p.jxy());
    double Jz = p.jz();
    double a = crop_01(force_scatter);
    if (alps::is_nonzero<1>(Jxy + std::abs(Jz))) {
      if (Jxy + Jz > 2 * a * Jxy) {
        // standard solutions
        v[1] = crop_0(std::min(Jxy/2, (Jz + Jxy)/4));
        v[2] = crop_0((Jz - Jxy)/2);
        v[3] = crop_0(std::min(Jxy/2, (-Jz + Jxy)/4));
        v[4] = crop_0(-(Jz + Jxy)/2);
      } else {
        // "ergodic" solutions
        v[1] = a*Jxy/2;
        v[2] = 0;
        v[3] = (1-a)*Jxy/2;
        v[4] = (-Jz + (2*a-1)*Jxy)/2;
      }
    } else {
      v[1] = v[2] = v[3] = v[4] = 0;
    }
    weight_ = std::max(v[1] + v[2], v[3] + v[4]);
    offset_ = Jz/4 - (v[1] + v[2]);
    sign_ = (p.jxy() >= 0 ? 1 : -1);
  }

  double density() const { return weight_; }
  double weight() const { return weight_; }
  double p_accept_para() const { return dip(v[1]+v[2], weight_); }
  double p_accept_anti() const { return dip(v[3]+v[4], weight_); }
  double p_accept(int r) const { return r ? p_accept_anti() : p_accept_para(); }
  double p_accept(int c0, int c1) const { return p_accept(c0 ^ c1); }
  double p_freeze_para() const { return dip(v[2], v[1]+v[2]); }
  double p_freeze_anti() const { return dip(v[4], v[3]+v[4]); }
  double p_freeze(int r) const { return r ? p_freeze_anti() : p_freeze_para(); }
  double p_freeze(int c0, int c1) const { return p_freeze(c0 ^ c1); }
  double p_reflect() const { return dip(v[3], v[1]+v[3]); }
  double offset() const { return offset_; }
  double sign() const { return sign_; }

  template<typename W>
  static bond_parameter check(const W& w)
  {
    using alps::is_equal;

    double rp = w.density() * w.p_accept_para();
    double w11 = rp * w.p_freeze_para();
    double w13 = rp - w11;
    double ra = w.density() * w.p_accept_anti();
    double w22 = ra * w.p_freeze_anti();
    double w23 = ra - w22;
    double w12 = -(w11+ w13 + w22 + w23) / 2;

    assert(w11 >= 0 && w13 >= 0 && w22 >= 0 && w23 >= 0);
    assert(alps::is_zero<1>(w23 - w.p_reflect() * (w13 + w23)));

    double jxy = 2 * (w13 + w23) * w.sign();
    double jz = 4 * (w12 + w11 + w13);

    assert(alps::is_zero<1>(w.offset() + w11 + w13 - jz/4));
    assert(alps::is_zero<1>(w.offset() + w22 + w23 - (-jz/4)));
    assert(alps::is_zero<1>(w.sign() * (w13 + w23) - jxy/2));

    return bond_parameter(0, jxy, jz);
  }

private:
  double weight_;
  double offset_;
  double sign_;
  double v[5]; // v[0] is not used
};


class site_weight {

  // loop equation and its solution:
  //
  //   w1 = Hx/2
  //
  //   density = w1

public:
  site_weight() : weight_(0), sign_(1) {}
  site_weight(const site_parameter& p)
  {
    weight_ = std::abs(p.hx()) / 2;
    sign_ = (p.hx() >= 0 ? 1 : -1);
  }

  void init(const site_parameter& p)
  {
    weight_ = std::abs(p.hx()) / 2;
    sign_ = (p.hx() >= 0 ? 1 : -1);
  }

  double density() const { return weight_; }
  double weight() const { return weight_; }
  double offset() const { return weight_; }
  double sign() const { return sign_; }

  template<typename W>
  static site_parameter check(const W& w)
  {
    double w1 = w.density();
    double hx = 2 * w1 * w.sign();
    assert(w1 >= 0);
    return site_parameter(0, hx);
  }

private:
  double density_;
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
