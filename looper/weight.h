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
#include <algorithm> // for std::min std::max
#include <cmath>     // for std::abs

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

namespace weight {

class xxz {

  // loop equations:
  //
  //   w12 + w11 + w13 = Jz/4
  //   w12 + w22 + w23 = -Jz/4
  //   w13 + w23 = |Jxy|/2
  //
  //   density = max((w11 + w13), (w22 + w23))
  //   Pa_para = (w11 + w13) / density
  //   Pa_anti = (w22 + w23) / density
  //   Pf_para = w11 / w13
  //   Pf_anti = w22 / w23
  //   Pr = w23 / (w13 + w 23)

  // standard solutions:
  //
  // i) Jz >= |Jxy|  (ferro-Ising)
  //      w12 = -Jz/4
  //      w11 = (Jz - |Jxy|)/2
  //      w13 = |Jxy|/2
  //      w22 = w23 = 0
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //      w12 = -|Jxy|/4
  //      w11 = w22 = 0
  //      w13 = (Jz + |Jxy|)/4
  //      w23 = (-Jz + |Jxy|)/4
  // iii) Jz <= -|Jxy|  (antiferro-Ising)
  //      w12 = Jz/4
  //      w11 = w13 = 0
  //      w22 = (-Jz - |Jxy|)/2
  //      w23 = |Jxy|/2

  // "ergodic" solutions (with additional parameter a)
  // i) Jz >= |Jxy|  (ferro-Ising)
  //      same as the standard solution
  // ii) -|Jxy| <= Jz <= |Jxy|  (XY)
  //   ii-1) Jz + |Jxy| >= 2 a |Jxy|
  //      same as the standard solution
  //   ii-2) Jz + |Jxy| <= 2 a |Jxy|
  //      w12 = (Jz - 2a |Jxy|)/4
  //      w11 = 0
  //      w22 = (-Jz + (2a-1) |Jxy|)/2
  //      w13 = a |Jxy| / 2
  //      w23 = (1-a) |Jxy| / 2
  // iii) Jz <= -|Jxy|  (antiferro-Ising)
  //      same as ii-2)

  //      w22 = (-Jz + a |Jxy|)/2
  //             (1-a) |Jxy| / 2

public:
  xxz() :
    density_(0), pa_para_(0), pa_anti_(0), pf_para_(0), pf_anti_(0), 
    p_reflect_(0), offset_(0), sign_(1) {}
  template<class P>
  xxz(const P& p, double force_scatter = 0) :
    density_(0), pa_para_(0), pa_anti_(0), pf_para_(0), pf_anti_(0),
    p_reflect_(0), offset_(0), sign_(1)
  {
    double Jxy = std::abs(p.jxy());
    double Jz = p.jz();
    double a = crop_01(force_scatter);
    double w12, w11, w13, w22, w23;
    if (!is_zero(Jxy + std::abs(Jz))) {
      if (Jxy + Jz > 2 * a * Jxy) {
        // standard solutions
        w12 = std::min(-std::abs(Jz)/4., -Jxy/4);
        w11 = std::max((Jz - Jxy)/2, 0.);
        w13 = std::max(std::min(Jxy/2, (Jz + Jxy)/4), 0.);
        w22 = std::max(-(Jz + Jxy)/2, 0.);
        w23 = std::max(std::min(Jxy/2, (-Jz + Jxy)/4), 0.);
      } else {
        // "ergodic" solutions
        w12 = (Jz - 2*a*Jxy)/4;
        w11 = 0;
        w13 = a*Jxy/2;
        w22 = (-Jz + (2*a-1)*Jxy)/2;
        w23 = (1-a)*Jxy/2;
      }

      density_ = std::max(w11+w13, w22+w23);
      pa_para_ = (w11+w13)/density_;
      pa_anti_ = (w22+w23)/density_;
      pf_para_ = (w11+w13>0) ? w11/(w11+w13) : 0;
      pf_anti_ = (w22+w23>0) ? w22/(w22+w23) : 0;
      p_reflect_ = (w13+w23>0) ? w23/(w13+w23) : 0;
      offset_ = Jz/4. - (w11+w13);
      sign_ = (p.jxy() >= 0 ? 1 : -1);
    }
  }

  double density() const { return density_; }
  double weight() const { return density_; }
  double p_accept_para() const { return pa_para_; }
  double p_accept_anti() const { return pa_anti_; }
  double p_accept(int r) const { return r ? pa_anti_ : pa_para_; }
  double p_accept(int c0, int c1) const { return p_accept(c0 ^ c1); }
  double p_freeze_para() const { return pf_para_; }
  double p_freeze_anti() const { return pf_anti_; }
  double p_freeze(int r) const { return r ? pf_anti_ : pf_para_; }
  double p_freeze(int c0, int c1) const { return p_freeze(c0 ^ c1); }
  double p_reflect() const { return p_reflect_; }
  double offset() const { return offset_; }
  double sign() const { return sign_; }

private:
  double density_;
  double pa_para_;
  double pa_anti_;
  double pf_para_;
  double pf_anti_;
  double p_reflect_;
  double offset_;
  double sign_;
};

template<typename P, typename W>
P check(const W& w)
{
  double rp = w.density() * w.p_accept_para();
  double w11 = rp * w.p_freeze_para();
  double w13 = rp - w11;
  double ra = w.density() * w.p_accept_anti();
  double w22 = ra * w.p_freeze_anti();
  double w23 = ra - w22;
  double w12 = -(w11+ w13 + w22 + w23) / 2;

  assert(w11 >= 0. && w13 >= 0. && w22 >= 0. && w23 >= 0.);
  assert(equal(w23, w.p_reflect() * (w13 + w23)));

  double jxy = 2 * (w13 + w23) * w.sign();
  double jz = 4 * (w12 + w11 + w13);

  // for sse
  assert(equal(w.offset() + w11 + w13, jz/4));
  assert(equal(w.offset() + w22 + w23, -jz/4));
  assert(equal(w.sign() * (w13 + w23), jxy/2));

  return P(0, jxy, jz);
}

} // end namespace weight


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
