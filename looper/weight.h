/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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
#include <algorithm> // for std::max
#include <cmath>     // for std::abs

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

namespace weight {

class ferro_ising
{
public:
  ferro_ising() : density_(0), offset_(0) {}
  template<class P>
  ferro_ising(const P& p) : density_(p.jz()/2), offset_(-p.jz()/4)
  {
    if (!(p.jx() == 0 && p.jy() == 0 && p.jz() >= 0))
      boost::throw_exception(std::invalid_argument(
        "invalid coupling constant"));
  }

  double density() const { return density_; }
  double weight() const { return density_; }
  static double p_accept_para() { return 1; }
  static double p_accept_anti() { return 0; }
  static double p_accept(int r) { return r ? 0 : 1; }
  static double p_accept(int c0, int c1) { return p_accept(c0 ^ c1); }
  static double p_freeze_para() { return 1; }
  static double p_freeze_anti() { return 0; }
#ifndef NDEBUG
  static double p_freeze(int r) {
    if (r != 0) boost::throw_exception(std::logic_error("ferro_ising"));
    return 1;
  }
  static double p_freeze(int c0, int c1) {
    if (c0 ^ c1 != 0) boost::throw_exception(std::logic_error("ferro_ising"));
    return 1;
  }
#else
  static double p_freeze(int) { return 1; }
  static double p_freeze(int, int) { return 1; }
#endif
  static double p_reflect() { return 0; }
  double offset() const { return offset_; }
  static double sign() { return 1; }

private:
  double density_;
  double offset_;
};

class antiferro_ising
{
public:
  antiferro_ising() : density_(0), offset_(0) {}
  template<class P>
  antiferro_ising(const P& p) : density_(-p.jz()/2), offset_(p.jz()/4)
  {
    if (!(p.jx() == 0 && p.jy() == 0 && p.jz() <= 0))
      boost::throw_exception(std::invalid_argument(
        "invalid coupling constant"));
  }

  double density() const { return density_; }
  double weight() const { return density_; }
  static double p_accept_para() { return 0; }
  static double p_accept_anti() { return 1; }
  static double p_accept(int r) { return r ? 1 : 0; }
  static double p_accept(int c0, int c1) { return p_accept(c0 ^ c1); }
  static double p_freeze_para() { return 0; }
  static double p_freeze_anti() { return 1; }
#ifndef NDEBUG
  static double p_freeze(int r) {
    if (r != 1) boost::throw_exception(std::logic_error("antiferro_ising"));
    return 1;
  }
  static double p_freeze(int c0, int c1) {
    if (c0 ^ c1 != 1)
      boost::throw_exception(std::logic_error("antiferro_ising"));
    return 1;
  }
#else
  static double p_freeze(int) { return 1; }
  static double p_freeze(int, int) { return 1; }
#endif
  static double p_reflect() { return 0; }
  double offset() { return offset_; }
  static double sign() { return 1; }

private:
  double density_;
  double offset_;
};

class xxz {
public:
  xxz() :
    density_(0), pa_para_(0), pa_anti_(0), pf_para_(0), pf_anti_(0), 
    p_reflect_(0), offset_(0), sign_(1) {}
  template<class P>
  xxz(const P& p) :
    density_(0), pa_para_(0), pa_anti_(0), pf_para_(0), pf_anti_(0),
    p_reflect_(0), offset_(0), sign_(1)
  {
    using std::abs; using std::max;
    double Jxy = abs(p.jxy());
    double Jz = p.jz();
    if (Jxy + abs(Jz) > 1.0e-10) {
      density_ = max(abs(Jz) / 2, (Jxy + abs(Jz)) / 4);
      pa_para_ = (Jxy + Jz) / (Jxy + abs(Jz));
      pa_anti_ = (Jxy - Jz) / (Jxy + abs(Jz));
      pf_para_ = pf_anti_ = 1 - Jxy / abs(Jz);
      p_reflect_ = (Jxy - Jz) / (2 * Jxy);
      offset_ = -max(abs(Jz) / 4, Jxy / 4);
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

template<typename W>
inline xxz_parameter check(const W& w)
{
  double jxy;
  double jz;
  return xxz_parameter(0, jxy, jz);
}

} // end namespace weight

typedef weight::xxz default_weight;


template<class W>
class uniform_bond_chooser
{
public:
  typedef W weight_type;

  uniform_bond_chooser() : n_(), weight_(), gw_() {}
  template<class G, class M>
  uniform_bond_chooser(const G& vg, const M& m) : n_(), weight_(), gw_()
  { init(vg, m); }
  template<class P>
  uniform_bond_chooser(const P& p) : n_(), weight_(), gw_() { init(p); }

  template<class G, class M>
  void init(const G& vg, const M& m)
  {
    assert(m.num_bond_types() == 1);
    n_ = double(boost::num_edges(vg.graph));
    weight_ = weight_type(m.uniform_bond());
    gw_ = n_ * weight_.weight();
  }
  template<class P>
  void init(const P& p) { init(p.virtual_graph, p.model); }

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

template<class W>
class bond_chooser
{
public:
  typedef W weight_type;

  bond_chooser() : weight_(), rc_(), gw_(0) {}
  template<class G, class M>
  bond_chooser(const G& vg, const M& m) : weight_(), rc_(), gw_(0)
  { init (vg, m); }
  template<class P>
  bond_chooser(const P& p) : weight_(), rc_(), gw_(0) { init(p); }

  template<class G, class M>
  void init(const G& vg, const M& m)
  {
    weight_.clear();
    gw_ = 0.0;
    if (boost::num_edges(vg.graph) > 0) {
      typename boost::graph_traits<typename G::graph_type>::edge_iterator
        ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(vg.graph);
           ei != ei_end; ++ei)
        weight_.push_back(
          weight_type(m.bond(boost::get(alps::edge_type_t(), vg.graph, *ei))));

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
  template<class P>
  void init(const P& p) { init(p.virtual_graph, p.model); }

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
