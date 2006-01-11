/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_GRAPH_H
#define LOOPER_GRAPH_H

#include "random_choice.h"
#include "union_find.h"

#include <alps/math.hpp>
#include <alps/osiris.h>
#include <boost/random.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <stdexcept>

namespace looper {

struct bond_graph_type {
  // g = 0 (g3 in textbook)
  //     1 (g4)
  //     2 (g1)
  //     3 (g2)
  static bool is_compatible(int g, int c0, int c1)
  { return (g >> 1) ^ c0 ^ c1; }
};

struct site_graph_type {
  // g = 0 (g1 in textbook)
  static bool is_compatible(int /* g */, int /* c */) { return true; }
};

template<class LOC>
class local_graph
{
public:
  typedef LOC location_t;

  local_graph() {}
  local_graph(int type, const location_t& loc) : type_(type), loc_(loc) {}

  int type() const { return type_; }
  bool is_compatible(int c) const
  { assert(is_site()); return site_graph_type::is_compatible(type_, c); }
  bool is_compatible(int c0, int c1) const
  { assert(is_bond()); return bond_graph_type::is_compatible(type_, c0, c1); }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  bool is_bond() const { return loc_.is_bond(); }
  bool is_site() const { return loc_.is_site(); }

  static local_graph bond_graph(int type, int pos)
  { return local_graph(type, location_t::bond_location(pos)); }
  static local_graph site_graph(int type, int pos)
  { return local_graph(type, location_t::site_location(pos)); }

  template<class T>
  boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr) const
  {
    assert(is_site());
    int loop1 = add(fragments);
    return boost::make_tuple(loop1, curr, loop1);
  }

  template<class T>
  boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */,
               int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr0, int curr1) const
  {
    assert(is_bond());
    int loop0, loop1;
    if (type_ == 0) {
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = add(fragments);
    } else if (type_ == 2) {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    } else {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    }
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }

  static int loop_l(int /* t */, int loop0, int /* loop1 */) { return loop0; }
  static int loop_u(int /* t */, int /* loop0 */, int loop1) { return loop1; }
  static int loop_l0(int /* t */, int loop0, int /* loop1 */) { return loop0; }
  static int loop_l1(int t, int loop0, int loop1) { return t ? loop1 : loop0; }
  static int loop_u0(int /* t */, int /* loop0 */, int loop1) { return loop1; }
  static int loop_u1(int t, int loop0, int loop1) { return t ? loop0 : loop1; }

private:
  int type_;
  location_t loc_;
};

template<class LOC>
inline int type(const local_graph<LOC>& g) { return g.type(); }
template<class LOC>
inline bool is_compatible(const local_graph<LOC>& g, int c)
{ return g.is_compatible(c); }
template<class LOC>
inline bool is_compatible(const local_graph<LOC>& g, int c0, int c1)
{ return g.is_compatible(c0, c1); }

template<class LOC>
inline int pos(const local_graph<LOC>& g) { return g.pos(); }
template<class LOC>
inline bool is_site(const local_graph<LOC>& g) { return g.is_site(); }
template<class LOC>
inline bool is_bond(const local_graph<LOC>& g) { return g.is_bond(); }

template<class T, class LOC>
boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph<LOC>& g, int curr)
{ return g.reconnect(fragments, curr); }

template<class T, class LOC>
boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph<LOC>& g,
          int curr0, int curr1)
{ return g.reconnect(fragments, curr0, curr1); }

// optimized version for antiferromagnetic Heisenberg interaction

template<class LOC>
class local_graph_haf
{
public:
  typedef LOC location_t;

  local_graph_haf() {}
  local_graph_haf(int type, const location_t& loc) : loc_(loc)
  { assert(type == 0); }

  static int type() { return 0; }
  bool is_compatible(int c) const
  {
    boost::throw_exception(std::logic_error("local_graph_haf::is_compatible"));
    return false;
  }
  bool is_compatible(int c0, int c1) const { return c0 ^ c1; }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  static bool is_bond() { return true; }
  static bool is_site() { return false; }

  static local_graph_haf bond_graph(int type, int pos)
  {
    assert(type == 0);
    return local_graph_haf(0, location_t::bond_location(pos));
  }
  static local_graph_haf site_graph(int type, int pos)
  {
    boost::throw_exception(std::logic_error("local_graph_haf::site_graph"));
    return local_graph_haf();
  }

  template<class T>
  boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& /* fragments */, int /* curr */) const
  {
    boost::throw_exception(std::logic_error("local_graph_haf::reconnect"));
    return boost::make_tuple(0, 0, 0);
  }

  template<class T>
  boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */,
               int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr0, int curr1) const
  {
    curr0 = unify(fragments, curr0, curr1);
    curr1 = add(fragments);
    return boost::make_tuple(curr1, curr1, curr0, curr1);
  }

  static int loop_l(int, int, int)
  {
    boost::throw_exception(std::logic_error("local_graph_haf::loop_u0"));
    return 0;
  }
  static int loop_u(int, int, int)
  {
    boost::throw_exception(std::logic_error("local_graph_haf::loop_u1"));
    return 0;
  }
  static int loop_l0(int, int loop0, int /* loop1 */) { return loop0; }
  static int loop_l1(int, int loop0, int /* loop1 */) { return loop0; }
  static int loop_u0(int, int /* loop0 */, int loop1) { return loop1; }
  static int loop_u1(int, int /* loop0 */, int loop1) { return loop1; }

private:
  location_t loc_;
};

template<class LOC>
inline int type(const local_graph_haf<LOC>& g) { return g.type(); }
template<class LOC>
inline bool is_compatible(const local_graph_haf<LOC>& g, int c)
{ return g.is_compatible(c); }
template<class LOC>
inline bool is_compatible(const local_graph_haf<LOC>& g, int c0, int c1)
{ return g.is_compatible(c0, c1); }

template<class LOC>
inline int pos(const local_graph_haf<LOC>& g) { return g.pos(); }
template<class LOC>
inline bool is_site(const local_graph_haf<LOC>& g) { return g.is_site(); }
template<class LOC>
inline bool is_bond(const local_graph_haf<LOC>& g) { return g.is_bond(); }

template<class T, class LOC>
boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph_haf<LOC>& g, int curr)
{ return g.reconnect(fragments, curr); }

template<class T, class LOC>
boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph_haf<LOC>& g,
          int curr0, int curr1)
{ return g.reconnect(fragments, curr0, curr1); }

//
// graph_chooser
//

template<class LOCAL_GRAPH, class ENGINE> class graph_chooser;

template<class LOC, class ENGINE>
class graph_chooser<local_graph<LOC>, ENGINE>
{
public:
  typedef local_graph<LOC>                   local_graph_t;
  typedef typename local_graph_t::location_t location_t;
  typedef ENGINE                             engine_t;

  graph_chooser(engine_t& eng)
    : r_uniform_(eng, boost::uniform_real<>(0, 1)),
      r_graph_(eng, random_choice<>()),
      r_time_(eng, boost::exponential_distribution<>()) {}
  template<class WEIGHT_TABLE>
  graph_chooser(engine_t& eng, const WEIGHT_TABLE& wt, double beta,
                bool is_path_integral)
    : r_uniform_(eng, boost::uniform_real<>(0, 1)),
      r_graph_(eng, random_choice<>()),
      r_time_(eng, boost::exponential_distribution<>())
  { init(wt, beta, is_path_integral); }

  template<class WEIGHT_TABLE>
  void init(const WEIGHT_TABLE& wt, double beta, bool is_path_integral)
  {
    graph_.resize(0);
    diag_.resize(0);
    offdiag_.resize(0);
    std::vector<double> w;
    weight_ = 0;
    for (typename WEIGHT_TABLE::site_weight_iterator
           itr = wt.site_weights().first;
         itr != wt.site_weights().second; ++itr) {
      for (int g = 0; g <= 0; ++g) {
        if (alps::is_nonzero<1>(itr->second.v[g])) {
          graph_.push_back(local_graph_t::site_graph(g, itr->first));
          w.push_back(itr->second.v[g]);
          weight_ += itr->second.v[g];
        }
      }
    }
    for (typename WEIGHT_TABLE::bond_weight_iterator
           itr = wt.bond_weights().first;
         itr != wt.bond_weights().second; ++itr) {
      for (int g = 0; g <= 3; ++g) {
        if (alps::is_nonzero<1>(itr->second.v[g])) {
          graph_.push_back(local_graph_t::bond_graph(g, itr->first));
          w.push_back(itr->second.v[g]);
          weight_ += itr->second.v[g];
        }
      }
      offdiag_.push_back(
        alps::is_nonzero<1>(itr->second.v[0] + itr->second.v[2]) ?
        itr->second.v[0] / (itr->second.v[0] + itr->second.v[2]) : 1);
    }
    if (!is_path_integral) {
      for (typename WEIGHT_TABLE::bond_weight_iterator
             itr = wt.bond_weights().first;
           itr != wt.bond_weights().second; ++itr) {
        diag_.push_back(
          alps::is_nonzero<1>(itr->second.v[0] + itr->second.v[1]) ?
            itr->second.v[0] / (itr->second.v[0] + itr->second.v[1]) : 1);
        diag_.push_back(
          alps::is_nonzero<1>(itr->second.v[2] + itr->second.v[3]) ?
            itr->second.v[2] / (itr->second.v[2] + itr->second.v[3]) : 1);
      }
    }
    if (w.size()) r_graph_.distribution().init(w);
    if (alps::is_zero<2>(weight_)) weight_ = 1.0e-20;
    if (is_path_integral)
      r_time_.distribution() =
	boost::exponential_distribution<>(weight_ / beta);
  }

  const local_graph_t& graph() const { return graph_[r_graph_()]; }
  local_graph_t diagonal(const location_t& loc, int /* c */) const
  { return local_graph_t(0, loc); }
  local_graph_t diagonal(const location_t& loc, int c0, int c1) const
  {
    int c = 1 ^ c0 ^ c1; // 0 for antiparallel, 1 for parallel
    int g = ((r_uniform_() < diag_[(pos(loc) << 1) | c]) ? 0 : 1) ^ (c << 1);
    return local_graph_t(g, loc);
  }
  local_graph_t offdiagonal(const location_t& loc) const
  {
    int g = (is_site(loc) || r_uniform_() < offdiag_[pos(loc)]) ? 0 : 2;
    return local_graph_t(g, loc);
  }
  double advance() const { return r_time_(); }
  double weight() const { return weight_; }

private:
  mutable boost::variate_generator<engine_t&,
    boost::uniform_real<> > r_uniform_;
  mutable boost::variate_generator<engine_t&,
    random_choice<> > r_graph_;
  mutable boost::variate_generator<engine_t&,
    boost::exponential_distribution<> > r_time_;
  double weight_;
  std::vector<local_graph_t> graph_;
  std::vector<double> diag_;
  std::vector<double> offdiag_;
};

template<class LOC, class ENGINE>
class graph_chooser<local_graph_haf<LOC>, ENGINE>
{
public:
  typedef local_graph_haf<LOC>                   local_graph_t;
  typedef typename local_graph_t::location_t     location_t;
  typedef ENGINE                                 engine_t;

  template<class WEIGHT_TABLE>
  graph_chooser(engine_t& eng, const WEIGHT_TABLE& wt, double beta,
                bool is_path_integral = true)
    : r_graph_(eng, random_choice<>()),
      r_time_(eng, boost::exponential_distribution<>())
  { init(wt, beta, is_path_integral); }

  template<class WEIGHT_TABLE>
  void init(const WEIGHT_TABLE& wt, double beta, bool is_path_integral)
  {
    diag_.claer();
    std::vector<double> w;
    weight_ = 0;
    for (typename WEIGHT_TABLE::bond_iterator itr = wt.bond_begin();
         itr != wt.bond_end(); ++itr) {
      if (alps::is_nonzero<1>(itr->second.v[0])) {
        diag_.push_back(local_graph_t::bond_graph(0, itr->first));
        w.push_back(itr->second.v[0]);
        weight_ += itr->second.v[0];
      }
    }
    if (w.size()) r_graph_.distribution().init(w);
    if (alps::is_zero<2>(weight_)) weight_ = 1.0e-20;
    if (is_path_integral)
      r_time_.distribution() =
        boost::exponential_distribution<>(weight_ / beta);
  }

  const local_graph_t& graph() const { return diag_[r_graph_()]; }
  static local_graph_t diagonal(const location_t& loc, int)
  { return local_graph_t(0, loc); }
  static local_graph_t diagonal(const location_t& loc, int, int)
  { return local_graph_t(0, loc); }
  static local_graph_t offdiagonal(const location_t& loc)
  { return local_graph_t(0, loc); }
  double advance() const { return r_time_(); }
  double weight() const { return weight_; }

private:
  mutable boost::variate_generator<engine_t&,
    random_choice<> > r_graph_;
  mutable boost::variate_generator<engine_t&,
    boost::exponential_distribution<> > r_time_;
  double weight_;
  std::vector<local_graph_t> diag_;
};

} // end namespace looper

#endif // LOOPER_GARPH_H
