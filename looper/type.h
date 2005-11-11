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

#include <alps/osiris.h>
#include <looper/union_find.h>
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
  location(int pos = 0, bool is_site = true) : loc_(pos << 1 + is_site) {}
  int pos() const { return loc_ >> 1; }
  bool is_site() const { return loc_ && 1; }
  bool is_bond() const { return !is_site(); }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
private:
  int loc_;
};

inline int pos(const location& loc) { return loc.pos(); }
inline bool is_site(const location& loc) { return loc.is_site(); }
inline bool is_bond(const location& loc) { return loc.is_bond(); }

inline location site_location(int pos) { return location(pos, false); }
inline location bond_location(int pos) { return location(pos, true); }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dump, const looper::location& loc)
{ loc.save(dump); return dump; }

inline alps::IDump& operator>>(alps::IDump& dump, looper::location& loc)
{ loc.load(dump); return dump; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

//
// graph
//

struct site_graph_type {
  /* g = 0 (g1 in textbook)
         1
         2 */
  static bool is_compatible(int g, int c) { return 2 - (c + g); }
  static bool is_locked(int g) { return g /* g != 0 */; }
};

struct bond_graph_type {
  /* g = 0 (g3 in textbook)
         1 (g4)
         2 (g1)
         3 (g2) */
  static bool is_compatible(int g, int c0, int c1)
  { return (g >> 1) ^ (c0 ^ c1); }
};

class local_graph
{
public:
  local_graph(int type, const location& loc) : type_(type), loc_(loc) {}

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

  const location& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  bool is_site() const { return loc_.is_site(); }
  bool is_bond() const { return loc_.is_bond(); }

  template<class T>
  boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr) const
  {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("is_compatible"));
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
    if (!is_site())
      boost::throw_exception(std::logic_error("is_compatible"));
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
  location loc_;
};

inline int type(const local_graph& g) { return g.type(); }
inline bool is_compatible(const local_graph& g, int c0, int c1)
{ return g.is_compatible(c0, c1); }
inline bool is_compatible(const local_graph& g, int c)
{ return g.is_compatible(c); }

inline int pos(const local_graph& g) { return g.pos(); }
inline bool is_site(const local_graph& g) { return g.is_site(); }
inline bool is_bond(const local_graph& g) { return g.is_bond(); }

inline local_graph site_graph(int g, int i)
{ return local_graph(g, site_location(i)); }
inline local_graph bond_graph(int g, int i)
{ return local_graph(g, bond_location(i)); }

template<class T>
boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph& g, int curr)
{ return g.reconnect(fragments, curr); }

template<class T>
boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph& g, int curr0, int curr1)
{ return g.reconnect(fragments, curr0, curr1); }

//
// local operator
//

struct local_operator_type
{
  BOOST_STATIC_CONSTANT(int, diagonal    = 0 /* 00 */);
  BOOST_STATIC_CONSTANT(int, offdiagonal = 1 /* 00 */);
  BOOST_STATIC_CONSTANT(int, identity    = 2 /* 10 */);
};

template<class QMC>
class local_operator;

template<>
class local_operator<sse> {
public:
  local_operator() {}
  local_operator(int type, const location& loc) : type_(type), loc_(loc) {}
  local_operator(const local_graph& g)
    : type_(local_operator_type::diagonal & (g.type() << 2)), loc_(g.loc()) {}

  void flip() { type_ ^= 1; }
  int type() const { return type_ & 3; }
  bool is_diagonal() const { return type() == local_operator_type::diagonal; }
  bool is_offdiagonal() const
  { return type() == local_operator_type::offdiagonal; }
  bool is_identity() const { return type() == local_operator_type::identity; }

  local_operator& operator=(const local_graph& g)
  {
    type_ = local_operator_type::diagonal & (g.type() << 2);
    loc_ = g.loc();
    return *this;
  }
  local_graph graph() const { return local_graph(graph_type(), loc_); }
  void assign_graph(const local_graph g) { type_ = type() & (g.type() << 2); }
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

  const location& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  bool is_site() const { return loc_.is_site(); }
  bool is_bond() const { return loc_.is_bond(); }

  void save(alps::ODump& dump) const { dump << type_ << loc_; }
  void load(alps::IDump& dump) { dump >> type_ >> loc_; }

private:
  int type_;
  location loc_;

public:
  int loop0, loop1; // no checkpoint
};

template<>
class local_operator<path_integral> : public local_operator<sse>
{
public:
  typedef local_operator<sse> super_type;
  local_operator() : super_type(), time_() {}
  local_operator(int type, const location& loc); // not defined
  local_operator(int type, const location& loc, double time)
    : super_type(type, loc), time_(time) {}
  local_operator(const local_graph& g); // not defined
  local_operator(const local_graph& g, double time)
    : super_type(g), time_(time) {}

  bool is_identity() const; // not defined

  double time() const { return time_; }

  void save(alps::ODump& dump) const { super_type::save(dump); dump << time_; }
  void load(alps::IDump& dump) { super_type::load(dump); dump >> time_; }

private:
  double time_;
};

template<class QMC>
int type(const local_operator<QMC>& op) { return op.type(); }
template<class QMC>
bool is_diagonal(const local_operator<QMC>& op) { return op.is_diagonal(); }
template<class QMC>
bool is_offdiagonal(const local_operator<QMC>& op)
{ return op.is_offdiagonal(); }
inline
bool is_identity(const local_operator<sse>& op) { return op.is_identity(); }

template<class QMC>
int pos(const local_operator<QMC>& op) { return op.pos(); }
template<class QMC>
bool is_site(const local_operator<QMC>& op) { return op.is_site(); }
template<class QMC>
bool is_bond(const local_operator<QMC>& op) { return op.is_bond(); }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class QMC>
alps::ODump& operator<<(alps::ODump& dump,
                        const looper::local_operator<QMC>& op)
{ op.save(dump); return dump; }

template<class QMC>
alps::IDump& operator>>(alps::IDump& dump, looper::local_operator<QMC>& op)
{ op.load(dump); return dump; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_NODE_H
