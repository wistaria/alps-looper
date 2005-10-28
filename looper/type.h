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

#ifndef LOOPER_TYPE_H
#define LOOPER_TYPE_H

#include <alps/osiris.h>

namespace looper {

//
// QMC type
//

struct path_integral {};
struct sse {};

enum operator_type { diagonal    = 0 /* 00 */, 
                     offdiagonal = 1 /* 01 */,
                     identity    = 2 /* 10 */};

template<class QMC> operator_graph_base;

template<>
class operator_graph_base<sse> {
public:
  operator_graph_base(operator_type type, unsigned int loc, bool is_bond)
    : type_(type), loc_(loc << 1 + is_bond) {}

  void flip() { type_ ^= 1; }
  operator_type type() const { return type_; }
  unsigned int location() const { return loc_ >> 1; }
  bool is_bond() const { return loc_ & 1; }
  bool is_site() const { return !is_bond(); }

  unsigned int loop_0() const { return loop0_; }
  unsigned int loop_1() const { return loop1_; }
  unsigned int set_loop_0(unsigned int loop) const { return loop0_ = loop; }
  unsigned int set_loop_1(unsigned int loop) const { return loop1_ = loop; }

  void save(alps::ODump& od) const { od << type_ << loc_; }
  void load(alps::IDump& id) { id >> type_ >> loc_; }

private:
  operator_type type_;
  unsigned int loc_;
  unsigned int loop0_, loop1_;
};

template<>
class operator_graph_base<path_integral> : public operator_graph_base<sse>
{
public:
  typedef operator_graph_base<sse> super_type;
  operator_graph(unsigned int type, unsigned int loc, bool is_bond,
		 double time)
    : super_type(type, loc, is_bond), time_(time) {}
  double time() const { return time_; }
  void save(alps::ODump& dp) const { super_type::save(dp); dp << time_; }
  void load(alps::IDump& dp) { super_type::load(dp); dp >> time_; }
private:
  double time_;
};


//
// graph
//

template<class QMC> 
class local_graph : public operator_graph_base<QMC>
{
public:
  typedef operator_graph_base<QMC> super_type;
  local_graph(unsigned_int type, unsigned int loc, bool is_bond)
    : super_type(type, loc, is_bond) {}
  
  static local_graph bond_graph(unsigned int type, unsigned int loc)
  {
#ifndef NDEBUG
    assert(type >= 1 && type <= 4);
#endif
    return local_graph(type, loc, true);
  }
  static local_graph site_graph(unsigned int type, unsigned int loc)
  {
#ifndef NDEBUG
  assert(type >= 1 && type <= 3);
#endif
  return local_graph(type, loc, false);
  }
};


// bond graph
inline
bool is_compatible(const local_graph& g, unsigned int c0, unsigned int c1)
{
#ifndef NDEBUG
  assert(g.is_bond());
#endif
  return !(((g.type()-1) >> 1) ^ (c0 ^ c1));
}

// site graph
inline
bool is_compatible(const local_graph& g, unsigned int c)
{
#ifndef NDEBUG
  assert(g.is_site());
#endif
  return (g.type() == 1) || ((3-g.type()) ^ c);
}

//
// operator
//

inline local_operator<sse> identity_operator()
{ return local_operator<sse>(true, identity, 0); }

inline local_operator<sse> site_diagonal_operator(unsigned int loc)
{ return local_operator<sse>(false, diagonal, loc); }

inline local_operator<sse> bond_diagonal_operator(unsigned int loc)
{ return local_operator<sse>(true, diagonal, loc); }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class QMC>
alps::ODump& operator<<(alps::ODump& dp, const looper::local_operator<QMC>& op)
{ op.save(dp); return dp; }

template<class QMC>
alps::IDump& operator>>(alps::IDump& dp, looper::local_operator<QMC>& op)
{ op.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_NODE_H
