/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: weight.h 554 2003-11-12 02:36:24Z wistaria $

#ifndef LOOPER_WEIGHT_H
#define LOOPER_WEIGHT_H

#include <looper/random_choice.h>

#include <algorithm> // for std::max, std::min
#include <cmath> // for std::abs

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

class default_weight
{
public:
  default_weight() : density_(0), p_freeze_(0), pa_para_(0),
		     pa_anti_(0), p_reflect_(0), offset_(0) {}
  template<class P>
  default_weight(const P& p) : density_(0), p_freeze_(0), pa_para_(0),
			       pa_anti_(0), p_reflect_(0), offset_(0)
  {
    using std::abs; using std::max;
    double Jxy = abs(p.jxy()); // ignore negative signs
    double Jz = p.jz();
    if (Jxy + abs(Jz) > 1.0e-10) {
      density_ = max(abs(Jz) / 2, (Jxy + abs(Jz)) / 4);
      p_freeze_ = range_01(1 - Jxy / abs(Jz));
      pa_para_ = range_01((Jxy + Jz) / (Jxy + abs(Jz)));
      pa_anti_ = range_01((Jxy - Jz) / (Jxy + abs(Jz)));
      p_reflect_ = range_01((Jxy - Jz) / (2 * Jxy));
      offset_ = -max(abs(Jz) / 4, Jxy / 4);
    }
  }

  double density() const { return density_; }
  double weight() const { return density_; }
  double p_freeze() const { return p_freeze_; }
  double p_accept_para() const { return pa_para_; }
  double p_accept_anti() const { return pa_anti_; }
  double p_accept(int c0, int c1) const {
    return (c0 ^ c1) ? pa_anti_ : pa_para_;
  }
  double p_reflect() const { return p_reflect_; }
  double offset() const { return offset_; }

protected:
  double range_01(double x) const {
    return std::min(std::max(x, 0.), 1.);
  }

private:
  double density_;
  double p_freeze_;
  double pa_para_;
  double pa_anti_;
  double p_reflect_;
  double offset_;
};
  

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

  weight_type& weight(int) { return w_; }
  const weight_type& weight(int) const { return w_; }
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
	  weight_type(m.bond(boost::get(edge_type_t(), vg.graph, *ei))));

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
