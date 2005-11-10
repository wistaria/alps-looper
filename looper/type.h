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
#include <looper/union_find.h>

namespace looper {

//
// QMC type
//

struct path_integral {};
struct sse {};

struct local_operator_type
{
  BOOST_STATIC_CONSTANT(unsigned int, diagonal    = 0 /* 00 */);
  BOOST_STATIC_CONSTANT(unsigned int, offdiagonal = 1 /* 00 */);
  BOOST_STATIC_CONSTANT(unsigned int, identity    = 2 /* 10 */);
};

template<class QMC>
class local_operator;

template<>
class local_operator<sse> {
public:
  local_operator() {}
  local_operator(unsigned int type, unsigned int loc, bool is_bond)
    : type_(type), loc_(loc << 1 + is_bond) {}

  void flip() { type_ ^= 1; }
  unsigned int type() const { return type_; }
  unsigned int location() const { return loc_ >> 1; }
  bool is_bond() const { return loc_ & 1; }
  bool is_site() const { return !is_bond(); }

  unsigned int loop_0() const { return loop0_; }
  unsigned int loop_1() const { return loop1_; }
  unsigned int set_loop_0(unsigned int loop) { return loop0_ = loop; }
  unsigned int set_loop_1(unsigned int loop) { return loop1_ = loop; }

  void save(alps::ODump& od) const { od << type_ << loc_; }
  void load(alps::IDump& id) { id >> type_ >> loc_; }

private:
  unsigned int type_;
  unsigned int loc_;
  unsigned int loop0_, loop1_;
};

template<>
class local_operator<path_integral> : public local_operator<sse>
{
public:
  typedef local_operator<sse> super_type;
  local_operator() : super_type() {}
  local_operator(unsigned int type, unsigned int loc, bool is_bond,
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

struct local_graph_type {
  BOOST_STATIC_CONSTANT(unsigned int, bond_g1 = 0);
  BOOST_STATIC_CONSTANT(unsigned int, bond_g2 = 1);
  BOOST_STATIC_CONSTANT(unsigned int, bond_g3 = 2);
  BOOST_STATIC_CONSTANT(unsigned int, bond_g4 = 3);
  BOOST_STATIC_CONSTANT(unsigned int, site_g1 = 4);
  BOOST_STATIC_CONSTANT(unsigned int, site_g2 = 5);
  BOOST_STATIC_CONSTANT(unsigned int, site_g3 = 6);

  static unsigned int bond_graph(unsigned int g) { return g-1; }
  static unsigned int site_graph(unsigned int g) { return g+3; }

  static bool is_bond(unsigned int g) { return !is_site(g); }
  static bool is_site(unsigned int g) { return g & 4; }
  static bool is_compatible(unsigned int g, unsigned int c0, unsigned int c1)
  {
#ifndef NDEBUG
    assert(is_bond(g));
#endif
    return !((g >> 1) ^ (c0 ^ c1));
  }
  static bool is_compatible(unsigned int g, unsigned int c)
  {
#ifndef NDEBUG
    assert(is_site(g));
#endif
    return ((g == site_g1) || ((g & 1) ^ c));
  }
};

class local_graph
{
public:
  local_graph(unsigned int type, unsigned int loc)
    : type_(type), loc_(loc) {}
  unsigned int type() const { return type_; }
  unsigned int location() const { return loc_; }
  bool is_bond() const { return local_graph_type::is_bond(type_); }
  bool is_site() const { return local_graph_type::is_site(type_); }
private:
  unsigned int type_, loc_;
};

inline
unsigned int bond_graph(unsigned int g)
{ return local_graph_type::bond_graph(g); }

inline
local_graph bond_graph(unsigned int g, unsigned int loc)
{ return local_graph(bond_graph(g), loc); }

inline
unsigned int site_graph(unsigned int g)
{ return local_graph_type::site_graph(g); }

inline
local_graph site_graph(unsigned int g, unsigned int loc)
{ return local_graph(site_graph(g), loc); }

inline
bool is_bond(unsigned int g) { return local_graph_type::is_bond(g); }

inline
bool is_bond(const local_graph& g) { return g.is_bond(); }

inline
bool is_site(unsigned int g) { return local_graph_type::is_site(g); }

inline
bool is_site(const local_graph& g) { return g.is_site(); }

inline
bool is_compatible(unsigned int g, unsigned int c0, unsigned int c1)
{ return local_graph_type::is_compatible(g, c0, c1); }

inline
bool is_compatible(const local_graph& g, unsigned int c0, unsigned int c1)
{ return local_graph_type::is_compatible(g.type(), c0, c1); }

inline
bool is_compatible(const local_graph& g, unsigned int c)
{ return local_graph_type::is_compatible(g.type(), c); }

inline
bool is_compatible(unsigned int g, unsigned int c)
{ return local_graph_type::is_compatible(g, c); }

//
// cluster fragment
//

typedef looper::union_find::node_idx clsuter_fragment;

//
// cluster_info
//

template<class QMC> struct cluster_info;

template<>
struct cluster_info<path_integral> {
  cluster_info(bool t = false)
    : to_flip(t), mag0(0), size(0), mag(0), length(0) {}
  bool to_flip;
  int mag0;
  int size;
  double mag;
  double length;
};

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
