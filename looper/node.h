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

#ifndef LOOPER_NODE_H
#define LOOPER_NODE_H

#include <looper/union_find.h>

#include <alps/osiris.h>
#include <boost/integer_traits.hpp>
// #include <boost/iterator.hpp>
#include <boost/static_assert.hpp>
#include <bitset>

namespace looper {

class local_graph {
public:
  local_graph(bool is_bond, unsigned int loc, unsigned int type)
    : loc_(is_bond ? 2*loc+1 : 2*loc), type_(type) {}

  unsigned int location() const { return loc_ >> 1; }
  bool is_bond() const { return loc_ & 1; }
  bool is_site() const { return !is_bond(); }
  unsigned int type() const { return type_; }

private:
  unsigned int loc_;
  unsigned int type_;
};

inline local_graph bond_graph(unsigned int loc, unsigned int type)
{
#ifndef NDEBUG
  assert(type >= 1 && type <= 4);
#endif
  return local_graph(true, loc, type);
}

inline local_graph site_graph(unsigned int loc, unsigned int type)
{
#ifndef NDEBUG
  assert(type >= 1 && type <= 3);
#endif
  return local_graph(false, loc, type);
}

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

namespace detail {

struct bits
{
  // number of bits
  BOOST_STATIC_CONSTANT(unsigned int, N = 7);

  BOOST_STATIC_CONSTANT(unsigned int, CONF = 0); // bit denoting spin direction

  BOOST_STATIC_CONSTANT(unsigned int, REFL = 1);
  BOOST_STATIC_CONSTANT(unsigned int, FREZ = 2); // only for easy-axis cases
  BOOST_STATIC_CONSTANT(unsigned int, ANTI = 3); // only for XYZ cases

  // for path integral
  BOOST_STATIC_CONSTANT(unsigned int, ADDD = 4); // set for newly-added node

  // for SSE
  BOOST_STATIC_CONSTANT(unsigned int, IDNT = 5); // identity operator
  BOOST_STATIC_CONSTANT(unsigned int, DIAG = 6); // diagonal operator

  // bit mask for clear()
  BOOST_STATIC_CONSTANT(unsigned int,
                        M_CLEAR = (1 << CONF) | (1 << IDNT) | (1 << DIAG));
};

} // end namespace detail

class node_property
{
public:
  typedef detail::bits bits;
  typedef std::bitset<bits::N> bitset;
  typedef bitset::reference reference;

  node_property() : prop_(0) {}
  void reset() { prop_.reset(); }
  void clear_graph() { prop_ &= bitset(bits::M_CLEAR); }

  reference conf() { return prop_[bits::CONF]; }
  bool conf() const { return prop_[bits::CONF]; }
  void flip_conf() { prop_.flip(bits::CONF); }

  bool is_refl() const { return prop_.test(bits::REFL); }
  bool is_frozen() const { return prop_.test(bits::FREZ); }
  bool is_anti() const { return prop_.test(bits::ANTI); }

  // for path integral
  bool is_new() const { return prop_.test(bits::ADDD); }
  bool is_old() const { return !is_new(); }
  void set_new(unsigned int is_refl, unsigned int is_frozen) {
    prop_.set(bits::REFL, is_refl).set(bits::FREZ, is_frozen)
      .set(bits::ADDD);
  }
  void set_new(unsigned int is_refl, unsigned int is_frozen,
               unsigned int is_anti) {
    set_new(is_refl, is_frozen);
    prop_.set(bits::ANTI, is_anti);
  }
  void set_old(unsigned int is_refl) { prop_.set(bits::REFL, is_refl); }
  void set_old(unsigned int is_refl, unsigned int is_frozen) {
    set_old(is_refl);
    prop_.set(bits::FREZ, is_frozen);
  }

  // for SSE
  bool is_identity() const { return prop_[bits::IDNT]; }
  bool is_diagonal() const { return prop_[bits::DIAG]; }
  bool is_offdiagonal() const { return !is_identity() && !is_diagonal(); }
  void set_to_identity() { prop_.set(bits::IDNT).reset(bits::DIAG); }
  void identity_to_diagonal() {
#ifndef NDEBUG
    assert(is_identity());
#endif
    prop_.reset(bits::IDNT).set(bits::DIAG);
  }
  void diagonal_to_identity() {
#ifndef NDEBUG
    assert(is_diagonal());
#endif
    prop_.reset(bits::DIAG).set(bits::IDNT);
  }
  void diagonal_to_offdiagonal() {
#ifndef NDEBUG
    assert(is_diagonal());
#endif
    prop_.reset(bits::DIAG);
  }
  void offdiagonal_to_diagonal() {
#ifndef NDEBUG
    assert(is_offdiagonal());
#endif
    prop_.set(bits::DIAG);
  }
  void flip_operator() {
#ifndef NDEBUG
    assert(!is_identity());
#endif
    prop_.flip(bits::DIAG);
  }

  std::ostream& output(std::ostream& os) const {
    return os << "node_property = " << prop_;
  }
  alps::ODump& save(alps::ODump& od) const {
    return od << uint32_t(prop_.to_ulong());
  }
  alps::IDump& load(alps::IDump& id) {
    prop_ = bitset(uint32_t(id));
    return id;
  }

private:
  bitset prop_;
};

class qmc_node : public node_property
{
public:
  typedef int segment_type;

  qmc_node() : node_property(), bond_(0), segment0_(), segment1_() {}

  uint32_t bond() const { return bond_; }
  void set_bond(uint32_t b) { bond_ = b; }

  segment_type& loop_segment(int i) {
    return (i == 0 ? segment0_ : segment1_);
  }
  segment_type loop_segment(int i) const {
    return (i == 0 ? segment0_ : segment1_);
  }

  void clear_graph() {
    node_property::clear_graph();
  }

  std::ostream& output(std::ostream& os) const {
    return node_property::output(os) << " bond = " << bond_;
  }
  alps::ODump& save(alps::ODump& od) const {
    return node_property::save(od) << bond_;
    // segment[01]_ are not saved
  }
  alps::IDump& load(alps::IDump& id) {
    return node_property::load(id) >> bond_;
    // segment[01]_ are not restored
  }

private:
  uint32_t bond_;
  segment_type segment0_;
  segment_type segment1_;
};


// // for path-integral

// template<class Itr>
// inline typename Itr::value_type::segment_type&
// segment_d(const Itr& itr)
// {
//   if (itr.at_boundary()) {
//     return itr->loop_segment(0);
//   } else {
//     if (itr->is_refl()) {
//       return itr->loop_segment(1);
//     } else {
//       return itr->loop_segment(itr.leg());
//     }
//   }
// }

// template<class Itr>
// inline typename Itr::value_type::segment_type&
// segment_u(const Itr& itr)
// {
//   if (itr.at_boundary()) {
//     return itr->loop_segment(0);
//   } else {
//     if (itr->is_refl()) {
//       return itr->loop_segment(0);
//     } else {
//       return itr->loop_segment(1-itr.leg());
//     }
//   }
// }

// // for sse

// template<class Itr>
// inline typename boost::iterator_value<Itr>::type::segment_type&
// segment_d(const Itr& itr, int leg)
// {
//   if (itr->is_refl()) {
//     return itr->loop_segment(1);
//   } else {
//     return itr->loop_segment(leg);
//   }
// }

// template<class Itr>
// inline typename boost::iterator_value<Itr>::type::segment_type&
// segment_u(const Itr& itr, int leg)
// {
//   if (itr->is_refl()) {
//     return itr->loop_segment(0);
//   } else {
//     return itr->loop_segment(1-leg);
//   }
// }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& od, const looper::qmc_node& n)
{ return n.save(od); }

inline alps::IDump& operator>>(alps::IDump& id, looper::qmc_node& n)
{ return n.load(id); }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_NODE_H
