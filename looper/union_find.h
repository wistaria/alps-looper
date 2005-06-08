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

// Weighted Union-Find Algorithm
// Reference:
//   D. Knuth,
//   `The Art of Computer Programming, Vol. 1, Fundamental Algorithms'
//   3rd edition (Addison Wesley, Reading, 1997) Sec 2.3.3.

#ifndef LOOPER_UNION_FIND_H
#define LOOPER_UNION_FIND_H

#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <vector>

namespace looper {

namespace union_find {

//
// index based implementation
//

struct node_idx_base {
  node_idx_base() : parent(-1) {}
  bool is_root() const { return parent < 0; }
  int weight() const { return -parent; }
  int add_weight(const node_idx_base& f) const { return parent += f.parent; }
  mutable int parent;  // negative for root fragment
};

struct node_idx : public node_idx_base {
  node_idx() : node_idx_base(), id() {}
  unsigned int id;
};

template<class T>
int add_node(std::vector<T>& v)
{
  v.push_back(T());
  return v.size() - 1; // return index of new node
}

template<class T>
inline int root_index(const std::vector<T>& v, int g)
{ while (!v[g].is_root()) g = v[g].parent; return g; }

template<class T>
inline void update_link(const std::vector<T>& v, int g, int r)
{ while (g != r) { int p = v[g].parent; v[g].parent = r; g = p; } }

template<class T>
inline const T& root(const std::vector<T>& v, int g)
{ int r = root_index(v, g); update_link(v, g, r); return v[r]; }

template<class T>
inline int cluster_id(const std::vector<T>& v, int g)
{ return root(v, g).id; }

template<class T>
inline int unify(std::vector<T>& v, int g0, int g1)
{
  int r0 = root_index(v, g0);
  int r1 = root_index(v, g1);
  if (r0 != r1) {
    if (v[r0].weight() < v[r1].weight()) std::swap(r0, r1);
    v[r0].add_weight(v[r1]);
    v[r1].parent = r0;
  }
  update_link(v, g0, r0);
  update_link(v, g1, r0);
  return r0; // return (new) root node
}

//
// pointer based implementation
//

namespace detail {

struct null_node {
  void reset() {}
  null_node& operator+=(const null_node&) { return *this; }
};

} // end namespace detail

//
// class template node<>
//
// Argument of unify() function should be type of node<>.
//
// In the base class of node<>, `reset()' and `operator+=()' member
// functions, as well as default and copy constructors, should be
// defined, so that quantities of child cluster can be summed into
// properly.

template<class BASE = detail::null_node, class W = unsigned int>
class node : public BASE
{
public:
  typedef BASE base_type;
  typedef W    weight_type;

  node() : base_type(), parent_(0), weight_(1) {}
  node(const node& n) : base_type(n), parent_(n.is_root() ? 0 : n.parent_),
                        weight_(n.weight_) {}

  bool is_root() const { return parent_ == 0; }
  node* root() {
    node* r = this;
    while (!r->is_root()) r = r->parent();
    node* n = this;
    while (n != r) {
      node* p = n->parent();
      n->set_parent(r);
      n = p;
    }
    return r;
  }
  const node* root() const {
    const node* r = this;
    while (!r->is_root()) r = r->parent();
    const node* n = this;
    while (n != r) {
      const node* p = n->parent();
      n->set_parent(r);
      n = p;
    }
    return r;
  }

  void set_parent(const node* p) const { parent_ = const_cast<node*>(p); }
  node* parent() { return parent_; }
  const node* parent() const { return parent_; }

  weight_type weight() const {
#ifndef NDEBUG
    if (!is_root())
      boost::throw_exception(std::logic_error("union_find::node::weight() : not a root node"));
#endif
    return weight_;
  }

  node& operator+=(node& c) {
    base_type::operator+=(c);
    weight_ += c.weight();
    c.set_parent(this);
    return *this;
  }

  void reset() {
    base_type::reset();
    parent_ = 0;
    weight_ = weight_type(1);
  }

private:
  // root node     : 0
  // non-root node : points my parent node
  mutable node* parent_;

  // root node     : weight of cluster (number of nodes in the cluster)
  // non-root node : meaningless
  weight_type weight_;
};


//
// function unify()
//

template<class T>
inline bool unify(node<T>& node0, node<T>& node1)
{
  // NOTE: if both nodes belong to the same root, then return false

  typedef node<T> node_type;
  node_type* root0 = node0.root();
  node_type* root1 = node1.root();
  if (root0 == root1) {
    // both node belong to the same tree
    return false;
  } else {
    // both node belong to different trees
    if (root0->weight() >= root1->weight()) {
      *root0 += *root1;
    } else {
      *root1 += *root0;
    }
    return true;
  }
}

} // end namespace union_find

} // end namespace looper

#endif // LOOPER_UNION_FIND_H
