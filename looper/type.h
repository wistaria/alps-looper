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

#ifndef LOOPER_TYPE_H
#define LOOPER_TYPE_H

#include "random_choice.h"
#include "union_find.h"

#include <alps/math.hpp>
#include <alps/osiris.h>
#include <boost/random.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <stdexcept>

namespace looper {

//
// QMC type
//

struct path_integral {};
struct sse {};

//
// location
//

class location
{
public:
  location(int pos = 0, bool is_bond = true)
    : loc_(pos << 1 | (is_bond ? 1 : 0)) {}
  int pos() const { return loc_ >> 1; }
  bool is_bond() const { return loc_ & 1; }
  bool is_site() const { return !is_bond(); }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
  static location bond_location(int pos) { return location(pos, true); }
  static location site_location(int pos) { return location(pos, false); }
private:
  int loc_;
};

inline int pos(const location& loc) { return loc.pos(); }
inline bool is_bond(const location& loc) { return loc.is_bond(); }
inline bool is_site(const location& loc) { return loc.is_site(); }

// optimized version for models with bond terms only

class location_bond
{
public:
  location_bond(int pos = 0, bool is_bond = true) : loc_(pos)
  {
    if (!is_bond)
      boost::throw_exception(std::invalid_argument("location_bond"));
  }
  int pos() const { return loc_; }
  static bool is_bond() { return true; }
  static bool is_site() { return false; }
  static location_bond bond_location(int pos) { return location_bond(pos); }
  static location_bond site_location(int)
  {
    boost::throw_exception(std::logic_error("location_bond"));
    return location_bond();
  }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
private:
  int loc_;
};

inline int pos(const location_bond& loc) { return loc.pos(); }
inline bool is_bond(const location_bond&) { return true; }
inline bool is_site(const location_bond&) { return false; }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dump, const looper::location& loc)
{ loc.save(dump); return dump; }

inline alps::IDump& operator>>(alps::IDump& dump, looper::location& loc)
{ loc.load(dump); return dump; }

inline alps::ODump& operator<<(alps::ODump& dump,
                               const looper::location_bond& loc)
{ loc.save(dump); return dump; }

inline alps::IDump& operator>>(alps::IDump& dump, looper::location_bond& loc)
{ loc.load(dump); return dump; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

//
// local_graph
//

struct bond_graph_type {
  // g = 0 (g3 in textbook)
  //     1 (g4)
  //     2 (g1)
  //     3 (g2)
  static bool is_compatible(int g, int c0, int c1)
  { return (g >> 1) ^ (c0 ^ c1); }
};

struct site_graph_type {
  // g = 0 (g1 in textbook)
  //     1 (for up spin locked with ghost spin)
  //     2 (for down spin locked with ghost spin)
  static bool is_compatible(int g, int c) { return 2 - (c + g); }
  static bool is_locked(int g) { return g /* g != 0 */; }
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
  {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("is_compatible"));
#endif
    return site_graph_type::is_compatible(type_, c);
  }
  bool is_compatible(int c0, int c1) const
  {
#ifndef NDEBUG
    if (!is_bond())
      boost::throw_exception(std::logic_error("is_compatible"));
#endif
    return bond_graph_type::is_compatible(type_, c0, c1);
  }
  bool is_locked() const
  {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("is_locked"));
#endif
    return site_graph_type::is_locked(type_);
  }

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
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("reconnect"));
#endif
    int loop0, loop1;
    if (type_ == 0) {
      loop0 = curr;
      loop1 = curr = add(fragments);
    } else {
      loop0 = loop1 = curr;
    }
    return boost::make_tuple(curr, loop0, loop1);
  }

  template<class T>
  boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */,
               int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr0, int curr1) const
  {
#ifndef NDEBUG
    if (!is_bond())
      boost::throw_exception(std::logic_error("reconnect"));
#endif
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
  {
#ifndef NDEBUG
    if (type != 0)
      boost::throw_exception(std::invalid_argument("local_graph_haf"));
#endif
  }

  static int type() { return 0; }
  bool is_compatible(int c) const
  {
#ifndef NDEBUG
    boost::throw_exception(std::logic_error("local_graph_haf::is_compatible"));
#endif
    return false;
  }
  bool is_compatible(int c0, int c1) const { return c0 ^ c1; }
  bool is_locked() const
  {
#ifndef NDEBUG
    boost::throw_exception(std::logic_error("local_graph_haf::is_locked"));
#endif
    return false;
  }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  static bool is_bond() { return true; }
  static bool is_site() { return false; }

  static local_graph_haf bond_graph(int type, int pos)
  {
#ifndef NDEBUG
    if (type != 0)
      boost::throw_exception(std::invalid_argument("local_graph_haf::bond_graph"));
#endif
    return local_graph_haf(0, location_t::bond_location(pos));
  }
  static local_graph_haf site_graph(int g, int i)
  {
    boost::throw_exception(std::logic_error("local_graph_haf::site_graph"));
    return local_graph_haf();
  }

  template<class T>
  boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr) const
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
  graph_chooser(engine_t& eng, const WEIGHT_TABLE& wt)
    : r_uniform_(eng, boost::uniform_real<>(0, 1)),
      r_graph_(eng, random_choice<>()),
      r_time_(eng, boost::exponential_distribution<>())
  { init(wt); }


  template<class WEIGHT_TABLE>
  void init(const WEIGHT_TABLE& wt)
  {
    diag_.resize(0);
    offdiag_.resize(0);
    std::vector<double> w;
    double r = 0;
    for (typename WEIGHT_TABLE::site_iterator itr = wt.site_begin();
         itr != wt.site_end(); ++itr) {
      for (int g = 0; g <= 2; ++g) {
        if (alps::is_nonzero<1>(itr->second.v[g])) {
          diag_.push_back(local_graph_t::site_graph(g, itr->first));
          w.push_back(itr->second.v[g]);
          r += itr->second.v[g];
        }
      }
    }
    for (typename WEIGHT_TABLE::bond_iterator itr = wt.bond_begin();
         itr != wt.bond_end(); ++itr) {
      for (int g = 0; g <= 3; ++g) {
        if (alps::is_nonzero<1>(itr->second.v[g])) {
          diag_.push_back(local_graph_t::bond_graph(g, itr->first));
          w.push_back(itr->second.v[g]);
          r += itr->second.v[g];
        }
      }
      offdiag_.push_back(
        alps::is_nonzero<1>(itr->second.v[0] + itr->second.v[2]) ?
        itr->second.v[0] / (itr->second.v[0] + itr->second.v[2]) : 1);
    }
    r_graph_.distribution().init(w);
    r_time_.distribution() = boost::exponential_distribution<>(r);
  }

  local_graph_t diagonal() const { return diag_[r_graph_()]; }
  local_graph_t offdiagonal(const location_t& loc) const
  {
    int g = (is_site(loc) || r_uniform_() < offdiag_[pos(loc)]) ? 0 : 2;
    return local_graph_t(g, loc);
  }
  double advance() const { return r_time_(); }

private:
  mutable boost::variate_generator<engine_t&,
    boost::uniform_real<> > r_uniform_;
  mutable boost::variate_generator<engine_t&,
    random_choice<> > r_graph_;
  mutable boost::variate_generator<engine_t&,
    boost::exponential_distribution<> > r_time_;
  std::vector<local_graph_t> diag_;
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
  graph_chooser(engine_t& eng, const WEIGHT_TABLE& wt)
    : r_graph_(eng, random_choice<>()),
      r_time_(eng, boost::exponential_distribution<>())
  { init(wt); }

  template<class WEIGHT_TABLE>
  void init(const WEIGHT_TABLE& wt)
  {
    diag_.claer();
    std::vector<double> w;
    double r = 0;
    for (typename WEIGHT_TABLE::bond_iterator itr = wt.bond_begin();
         itr != wt.bond_end(); ++itr) {
      if (alps::is_nonzero<1>(itr->second.v[0])) {
        diag_.push_back(local_graph_t::bond_graph(0, itr->first));
        w.push_back(itr->second.v[0]);
        r += itr->second.v[0];
      }
    }
    r_graph_.distribution().init(w);
    r_time_.distribution() = boost::exponential_distribution<>(r);
  }

  local_graph_t diagonal() const { return diag_[r_graph_()]; }
  static local_graph_t offdiagonal(const location_t& loc)
  { return local_graph_t(0, loc); }
  double advance() const { return r_time_(); }

private:
  mutable boost::variate_generator<engine_t&,
    random_choice<> > r_graph_;
  mutable boost::variate_generator<engine_t&,
    boost::exponential_distribution<> > r_time_;
  std::vector<local_graph_t> diag_;
};

//
// local_operator
//

struct local_operator_type
{
  BOOST_STATIC_CONSTANT(int, diagonal    = 0 /* 00 */);
  BOOST_STATIC_CONSTANT(int, offdiagonal = 1 /* 00 */);
  BOOST_STATIC_CONSTANT(int, identity    = 2 /* 10 */);
};

template<class QMC, class LOC_G>
class local_operator;

template<class LOC_G>
class local_operator<sse, LOC_G>
{
public:
  typedef LOC_G                              local_graph_t;
  typedef typename local_graph_t::location_t location_t;

  local_operator() : type_(local_operator_type::identity), loc_() {}
  local_operator(int type, const location_t& loc) : type_(type), loc_(loc) {}
  local_operator(const local_graph_t& g)
    : type_(local_operator_type::diagonal & (g.type() << 2)), loc_(g.loc()) {}

  void flip() { type_ ^= 1; }
  int type() const { return type_ & 3; }
  bool is_diagonal() const { return type() == local_operator_type::diagonal; }
  bool is_offdiagonal() const
  { return type() == local_operator_type::offdiagonal; }
  bool is_identity() const { return type() == local_operator_type::identity; }

  local_graph_t graph() const { return local_graph_t(graph_type(), loc_); }
  void assign_graph(const local_graph_t g) { type_ = type() & (g.type() << 2); }
  void clear_graph() { type_ &= 3; }
  int graph_type() const { return type_ >> 2; }
  bool is_locked() const
  {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("is_locked"));
#endif
    return site_graph_type::is_locked(graph_type());
  }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  bool is_site() const { return loc_.is_site(); }
  bool is_bond() const { return loc_.is_bond(); }

  void save(alps::ODump& dump) const { dump << type_ << loc_; }
  void load(alps::IDump& dump) { dump >> type_ >> loc_; }

private:
  int type_;
  location_t loc_;

public:
  int loop0, loop1; // no checkpoint
};

template<class LOC_G>
class local_operator<path_integral, LOC_G> : public local_operator<sse, LOC_G>
{
public:
  typedef local_operator<sse, LOC_G>         super_type;
  typedef typename super_type::local_graph_t local_graph_t;
  typedef typename super_type::location_t    location_t;

  local_operator() : super_type(), time_() {}
  local_operator(int type, const location_t& loc); // not defined
  local_operator(int type, const location_t& loc, double time)
    : super_type(type, loc), time_(time) {}
  local_operator(const local_graph_t& g); // not defined
  local_operator(const local_graph_t& g, double time)
    : super_type(g), time_(time) {}

  bool is_identity() const; // not defined

  double time() const { return time_; }

  void save(alps::ODump& dump) const { super_type::save(dump); dump << time_; }
  void load(alps::IDump& dump) { super_type::load(dump); dump >> time_; }

private:
  double time_;
};

template<class QMC, class LOC_G>
int type(const local_operator<QMC, LOC_G>& op) { return op.type(); }
template<class QMC, class LOC_G>
bool is_diagonal(const local_operator<QMC, LOC_G>& op)
{ return op.is_diagonal(); }
template<class QMC, class LOC_G>
bool is_offdiagonal(const local_operator<QMC, LOC_G>& op)
{ return op.is_offdiagonal(); }
template<class LOC_G>
bool is_identity(const local_operator<sse, LOC_G>& op)
{ return op.is_identity(); }

template<class QMC, class LOC_G>
int pos(const local_operator<QMC, LOC_G>& op) { return op.pos(); }
template<class QMC, class LOC_G>
bool is_site(const local_operator<QMC, LOC_G>& op) { return op.is_site(); }
template<class QMC, class LOC_G>
bool is_bond(const local_operator<QMC, LOC_G>& op) { return op.is_bond(); }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class QMC, class LOC_G>
alps::ODump& operator<<(alps::ODump& dump,
                        const looper::local_operator<QMC, LOC_G>& op)
{ op.save(dump); return dump; }

template<class QMC, class LOC_G>
alps::IDump& operator>>(alps::IDump& dump,
                        looper::local_operator<QMC, LOC_G>& op)
{ op.load(dump); return dump; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_TYPE_H
