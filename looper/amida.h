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

#ifndef LOOPER_AMIDA_H
#define LOOPER_AMIDA_H

#include <alps/osiris.h>
#include <boost/throw_exception.hpp>
#include <boost/utility.hpp>
#include <deque>
#include <iterator>
#include <stdexcept>
#include <utility>

#include <looper/vector_helper.h>


//
// forward declarations
//

namespace looper {

struct amida_node_base
{
  // Convention
  //               link  cut bottom top   vacant
  //  series[0] :   s1    s1   s1    s1    max
  //  series[1] :   s2    s1   max   max   0
  //  next[0]   :   X1    X1   X1    0     X1       goal flag
  //  prev[0]   :   X2    X2   0     X1    X1       root flag
  //  next[1]   :   X3    X1   0     0     0        normal flag
  //  prev[1]   :   X4    0    0     0     0        link flag

  bool at_bottom() const { return prev[0] == 0; }
  bool at_top() const { return next[0] == 0; }
  bool at_boundary() const { return series[1] == _node_max; }
  bool is_link() const { return prev[1] != 0; }
  bool is_cut() const { return series[0] == series[1]; }
  bool is_node() const { return series[0] != _node_max; }
  bool is_vacant() const { return series[0] == _node_max; }

  void set_as_vacant(amida_node_base* t) {
    series[0] = _node_max;
    series[1] = 0;
    next[0] = t;
    prev[0] = t;
    next[1] = 0;
    prev[1] = 0;
  }

  void set_as_bottom(amida_node_base* t, std::size_t s) {
    series[0] = s;
    series[1] = _node_max;
    next[0] = t;
    prev[0] = 0;
    next[1] = 0;
    prev[1] = 0;
  }

  void set_as_top(amida_node_base* b, std::size_t s) {
    series[0] = s;
    series[1] = _node_max;
    next[0] = 0;
    prev[0] = b;
    next[1] = 0;
    prev[1] = 0;
  }

  std::size_t series[2];
  amida_node_base* next[2];
  amida_node_base* prev[2];

private:
  static const std::size_t _node_max = ~0; // = 111...111
};


template<class T>
struct amida_node : public amida_node_base
{
  amida_node() : amida_node_base(), data_() {}

  T data_;
  template<class H>
  void save(alps::ODump& od, const H& helper) const
  {
    od << series[0] << series[1]
       << ((next[0] == 0) ?
           0 : helper.index(static_cast<amida_node<T> *>(next[0])) + 1)
       << ((next[1] == 0) ?
           0 : helper.index(static_cast<amida_node<T> *>(next[1])) + 1)
       << ((prev[0] == 0) ?
           0 : helper.index(static_cast<amida_node<T> *>(prev[0])) + 1)
       << ((prev[1] == 0) ?
           0 : helper.index(static_cast<amida_node<T> *>(prev[1])) + 1)
       << data_;
  }
  template<class D>
  void load(alps::IDump& id, D& array)
  {
    std::size_t n0, n1, p0, p1;
    id >> series[0] >> series[1] >> n0 >> n1 >> p0 >> p1;
    next[0] = (n0 == 0 ? 0 : static_cast<amida_node_base *>(&array[n0 - 1]));
    next[1] = (n1 == 0 ? 0 : static_cast<amida_node_base *>(&array[n1 - 1]));
    prev[0] = (p0 == 0 ? 0 : static_cast<amida_node_base *>(&array[p0 - 1]));
    prev[1] = (p1 == 0 ? 0 : static_cast<amida_node_base *>(&array[p1 - 1]));
    id >> data_;
  }
};


struct amida_series_iterator_base
{
  typedef std::size_t                     size_type;
  typedef std::ptrdiff_t                  difference_type;
  typedef std::bidirectional_iterator_tag iterator_category;

  amida_node_base* node_;
  size_type ser_;
  size_type leg_;

  amida_series_iterator_base() : node_(), ser_(), leg_() {}
  amida_series_iterator_base(amida_node_base* x, size_type s)
    : node_(x), ser_(s), leg_()
  {
    if (node_)
      if (node_->series[0] == ser_)
        leg_ = 0;
      else
        leg_ = 1;
  }
  ~amida_series_iterator_base() {}

  void jump() { leg_ ^= 1; ser_ = node_->series[leg_]; }

  size_type series() const { return ser_; }
  size_type leg() const { return leg_; }

  bool at_boundary() const { return node_->at_boundary(); }
  bool at_bottom() const { return node_->at_bottom(); }
  bool at_top() const { return node_->at_top(); }

  void incr() {
    node_ = node_->next[leg_];
    if (node_->series[0] == ser_)
      leg_ = 0;
    else
      leg_ = 1;
  }
  void decr() {
    node_ = node_->prev[leg_];
    if (node_->series[0] == ser_)
      leg_ = 0;
    else
      leg_ = 1;
  }

  bool operator==(const amida_series_iterator_base& x) const
  { return node_ == x.node_ && ser_ == x.ser_; }

  bool operator!=(const amida_series_iterator_base& x) const
  { return node_ != x.node_ || ser_ != x.ser_; }
};


template<class T>
struct amida_const_series_iterator : public amida_series_iterator_base
{
  typedef std::size_t   size_type;
  typedef T             value_type;
  typedef const T&      const_reference;
  typedef const T*      const_pointer;
  typedef amida_node<T> node_type;
  typedef amida_const_series_iterator<T> self_;

  amida_const_series_iterator() {}
  amida_const_series_iterator(node_type* x, size_type s)
    : amida_series_iterator_base(x, s) {}

  const_reference operator*() const { return ((node_type*)node_)->data_; }
  const_pointer operator->() const { return &(operator*()); }

  self_ operator++() { this->incr(); return *this; }
  self_ operator++(int) { self_ tmp = *this; this->incr(); return tmp; }
  self_ operator--() { this->decr(); return *this; }
  self_ operator--(int) { self_ tmp = *this; this->decr(); return tmp; }
};


template<class T>
struct amida_series_iterator : public amida_const_series_iterator<T>
{
  typedef amida_const_series_iterator<T> base_;
  typedef typename base_::size_type      size_type;
  typedef T                              value_type;
  typedef T&                             reference;
  typedef T*                             pointer;
  typedef amida_node<T>                  node_type;
  typedef amida_series_iterator<T>       self_;

  amida_series_iterator() {}
  amida_series_iterator(node_type* x, size_type s) : base_(x, s) {}

  reference operator*() const { return ((node_type*)base_::node_)->data_; }
  pointer operator->() const { return &(operator*()); }

  self_ operator++() { this->incr(); return *this; }
  self_ operator++(int) { self_ tmp = *this; this->incr(); return tmp; }
  self_ operator--() { this->decr(); return *this; }
  self_ operator--(int) { self_ tmp = *this; this->decr(); return tmp; }
};


template<class T>
class amida
{
public:
  typedef std::size_t                    size_type;
  typedef T                              value_type;
  typedef T&                             reference;
  typedef const T&                       const_reference;
  typedef T*                             pointer;
  typedef const T*                       const_pointer;
  typedef amida_series_iterator<T>       iterator;
  typedef amida_const_series_iterator<T> const_iterator;
  typedef amida_node<T>                  node_type;

  // constructors
  explicit amida(size_type r = 0)
    : array_(), vacant_(), num_series_(), num_nodes_(), num_links_(),
    num_cuts_(), largest_(0) { init(r); }
  amida(const amida& a); // not implemented!
  ~amida() {}

  void init(size_type r)
  {
    num_series_ = r;
    array_.resize(2 * num_series_);
    for (size_type i = 0; i < num_series_; ++i) {
      array_[i              ].set_as_bottom(&array_[i + num_series_], i);
      array_[i + num_series_].set_as_top   (&array_[i              ], i);
    }

    // setup stack of vacant nodes
    largest_ = array_.size();
    vacant_ = 0;

    // set current number of nodes, etc.
    num_nodes_ = 2 * num_series_;
    num_links_ = 0;
    num_cuts_ = 0;
  }

  void clear() { init(num_series_); }

  size_type num_series() const { return num_series_; }
  size_type num_nodes() const { return num_nodes_; }
  size_type num_links() const { return num_links_; }
  size_type num_cuts() const { return num_cuts_; }
  size_type num_max_nodes() const { return largest_; }
  size_type capacity() const { return array_.capacity(); }
  double memory() const { return sizeof(node_type) * capacity(); }

  std::pair<iterator, iterator> series(size_type s)
  {
    return std::make_pair(iterator(&array_[s], s),
                          iterator(&array_[s + num_series_], s));
  }
  std::pair<const_iterator, const_iterator> series(size_type s) const
  {
    return std::make_pair(
      const_iterator(const_cast<node_type *>(&array_[s]), s),
      const_iterator(const_cast<node_type *>(&array_[s + num_series_]), s));
  }

  std::pair<iterator, iterator>
  insert_link(const T& t,
              const iterator& prev0, const iterator& prev1,
              const iterator& next0, const iterator& next1)
  {
    // get a new node
    node_type* n = static_cast<node_type *>(new_node());
    n->data_ = t;

    n->series[0] = prev0.ser_;
    n->series[1] = prev1.ser_;

    // adjust link pointers
    n->prev[0] = prev0.node_;
    n->prev[1] = prev1.node_;
    prev0.node_->next[prev0.leg_] = n;
    prev1.node_->next[prev1.leg_] = n;

    n->next[0] = next0.node_;
    n->next[1] = next1.node_;
    next0.node_->prev[next0.leg_] = n;
    next1.node_->prev[next1.leg_] = n;

    ++num_links_;
    return std::make_pair(iterator(n, n->series[0]),
                          iterator(n, n->series[1]));
  }

  std::pair<iterator, iterator>
  insert_link_next(const T& t, const iterator& prev0, const iterator& prev1) {
    return insert_link(t, prev0, prev1,
                       boost::next(prev0), boost::next(prev1));
  }

  std::pair<iterator, iterator>
  insert_link_prev(const T& t, const iterator& next0, const iterator& next1) {
    return insert_link(t, boost::prior(next0), boost::prior(next1),
                       next0, next1);
  }

  void erase(iterator itr)
  {
    if (itr.node_->is_link()) {
      iterator prev0, next0, prev1, next1;
      if (itr.ser_ == 0) {
        prev0 = boost::prior(itr);
        next0 = boost::next(itr);
        itr.jump();
        prev1 = boost::prior(itr);
        next1 = boost::next(itr);
      } else {
        prev1 = boost::prior(itr);
        next1 = boost::next(itr);
        itr.jump();
        prev0 = boost::prior(itr);
        next0 = boost::next(itr);
      }
      prev0.node_->next[prev0.leg_] = next0.node_;
      prev1.node_->next[prev1.leg_] = next1.node_;
      next0.node_->prev[next0.leg_] = prev0.node_;
      next1.node_->prev[next1.leg_] = prev1.node_;
      --num_links_;
    } else if (itr.node_->is_cut()) {
      iterator prev = boost::prior(itr);
      iterator next = boost::next(itr);
      prev.node_->next[prev.leg_] = next.node_;
      next.node_->prev[next.leg_] = prev.node_;
      --num_cuts_;
    } else {
      boost::throw_exception(std::logic_error("amida<>::erase() : couldn't erase node"));
    }

    delete_node(itr.node_);
  }

  void save(alps::ODump& od) const
  {
    index_helper<std::deque<node_type> > helper(array_);
    od << largest_ << num_series_ << num_nodes_ << num_links_ << num_cuts_
       << ((vacant_ == 0) ?
           0 : helper.index(static_cast<node_type *>(vacant_)) + 1);
    for (std::size_t i = 0; i < largest_; ++i) array_[i].save(od, helper);
  }
  void load(alps::IDump& id)
  {
    std::size_t vacant;
    id >> largest_ >> num_series_ >> num_nodes_ >> num_links_ >> num_cuts_
       >> vacant;
    array_.resize(largest_);
    vacant_ = static_cast<node_type*>(vacant == 0 ? 0 : &(array_[vacant - 1]));
    for (std::size_t i = 0; i < largest_; ++i) array_[i].load(id, array_);
  }

  // only for special purposes
  typedef std::deque<node_type> base_type;
  const std::deque<node_type>& base() const { return array_; }
  const amida_node_base * vacant() const { return vacant_; }

protected:
  amida_node_base * new_node()
  {
    ++num_nodes_;
    if (vacant_ == 0) {
      // stack expansion
      array_.resize(array_.size() + 1);
      amida_node_base * tmp = &(array_[array_.size() - 1]);
      ++largest_;
      return tmp;
    } else {
      // pop from stack
      amida_node_base * tmp = vacant_;
      vacant_ = vacant_->next[0];
      return tmp;
    }
  }

  void delete_node(amida_node_base * n)
  {
    // push to stack of vacant nodes
    n->set_as_vacant(vacant_);
    vacant_ = n;
    --num_nodes_;
  }

private:
  std::deque<node_type> array_;
  amida_node_base * vacant_;    // pointer to the top of vacant node stack
  size_type num_series_; // number of series
  size_type num_nodes_;  // number of nodes
  size_type num_links_;  // number of links
  size_type num_cuts_;   // number of cuts
  size_type largest_;    // max number of nodes
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class T>
alps::ODump& operator<<(alps::ODump& od, const looper::amida<T>& c) {
  c.save(od); return od;
}

template<class T>
alps::IDump& operator>>(alps::IDump& id, looper::amida<T>& c) {
  c.load(id); return id;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif


//
// only for debugging
//

namespace looper {

template<class T>
int index(const amida<T>& a, const typename amida<T>::iterator& itr)
{
  return index_helper<>::index(a.base(),
    static_cast<amida_node<T> *>(itr.node_));
}

template<class T>
void output_serial(const amida<T>& a, std::ostream& os = std::cout)
{
  os << "[amida Information]\n";
  for (std::size_t i = 0; i != a.base().size(); ++i) {
    amida_node<T> * p =
      const_cast<amida_node<T> *>(&(a.base()[i]));
    typename amida<T>::iterator s0(p, p->series[0]);
    typename amida<T>::iterator s1(p, p->series[1]);
    if (!p->is_vacant()) {
      if (p->at_boundary()) {
        if (p->at_bottom()) {
          os << index(a, s0) << '\t' << "root at " << s0.ser_;
          os << " :\t" << "next node is "
             << index(a, boost::next(s0));
          os << "  ";
          os << "node: " << index(a, s0) << ", contents: " << p->data_;
          os << std::endl;
        } else {
          os << index(a, s0) << '\t' << "goal at " << s0.ser_;
          os << " :\t" << "prev node is "
             << index(a, boost::prior(s0));
          os << "  ";
          os << "node: " << index(a, s0) << ", contents: " << p->data_;
          os << std::endl;
        }
      } else {
        os << index(a, s0) << '\t' << "node connecting "
           << s0.ser_ << " and " << s1.ser_;
        os << " :\t" << index(a, boost::next(s0));
        os << ' ' << index(a, boost::prior(s0));
        os << ' ' << index(a, boost::next(s1));
        os << ' ' << index(a, boost::prior(s1));
        os << '\t';
        os << "node: " << index(a, s0) << ", contents: " << p->data_;
        os << std::endl;
      }
    }
  }
}

template<class T>
void output_stack(const amida<T>& a, std::ostream& os = std::cout)
{
  os << "[amida Stack Information]\n";
  for (amida_node_base * s =
         const_cast<amida_node_base *>(a.vacant()); s != 0;
       s = s->next[0]) {
    typename amida<T>::iterator
      itr(static_cast<typename amida<T>::node_type *>(s),
          s->series[0]);
    os << index(a, itr) << "->";
  }
  os << "NULL\n";
}

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class T>
std::ostream& operator<<(std::ostream& os, const looper::amida<T>& a)
{ looper::output_serial(a, os); looper::output_stack(a, os); return os; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace looper
#endif

#endif // LOOPER_AMIDA_H
