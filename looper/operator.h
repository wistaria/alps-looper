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

#ifndef LOOPER_OPERATOR_H
#define LOOPER_OPERATOR_H

#include "type.h"

#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

struct local_operator_type
{
  BOOST_STATIC_CONSTANT(int, diagonal    = 0 /* 00 */);
  BOOST_STATIC_CONSTANT(int, offdiagonal = 1 /* 01 */);
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
  void assign_graph(const local_graph_t g)
  {
#ifndef NDEBUG
    if (!is_identity() && (loc_ != g.log()))
      boost::throw_exception(std::logic_error("assign_graph"));
#endif
    // if type==identity then type will be set to diagonal, otherwize unchanged
    type_ = (type_ & 1) | (g.type() << 2);
    loc_ = g.loc();
  }
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
