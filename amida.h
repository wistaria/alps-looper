/**************************************************************************** 
*
* alps/looper: multi-cluster quantum Monte Carlo algorithm for spin systems
*              in path-integral and SSE representations
*
* $Id: amida.h 408 2003-10-10 09:34:54Z wistaria $
*
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>,
*
* Permission is hereby granted, free of charge, to any person or organization 
* obtaining a copy of the software covered by this license (the "Software") 
* to use, reproduce, display, distribute, execute, and transmit the Software, 
* and to prepare derivative works of the Software, and to permit others
* to do so for non-commerical academic use, all subject to the following:
*
* The copyright notice in the Software and this entire statement, including 
* the above license grant, this restriction and the following disclaimer, 
* must be included in all copies of the Software, in whole or in part, and 
* all derivative works of the Software, unless such copies or derivative 
* works are solely in the form of machine-executable object code generated by 
* a source language processor.
*
* In any scientific publication based in part or wholly on the Software, the
* use of the Software has to be acknowledged and the publications quoted
* on the web page http://www.alps.org/license/ have to be referenced.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
**************************************************************************/

#ifndef LOOPER_AMIDA_H
#define LOOPER_AMIDA_H

#include <alps/osiris.h>
#include <boost/throw_exception.hpp>
#include <boost/type_traits.hpp> // for boost::is_class
#include <cmath>
#include <deque>
#include <iosfwd>
#include <iterator>
#include <limits>
#include <stdexcept>

//
// forward declarations
//

namespace looper {

template<class T, std::size_t C = 10> class amida;

namespace detail {

struct amida_node_base;
template<class T, bool IsClass> class amida_node;
template<class T> class amida_pointer;

} // end namespace detail

template<class T> class amida_series_iterator;

template<class T> std::size_t series(const amida_series_iterator<T>&);
template<class T> std::size_t leg(const amida_series_iterator<T>&);
template<class T> void jump(amida_series_iterator<T>&);
template<class T> void proceed(amida_series_iterator<T>&, size_t);
  // template<class T, std::size_t C> void rotate_series(amida<T, C>&, std::size_t, std::size_t, int);

namespace detail {

template<class T> inline amida_node_base* addr(const amida_pointer<T>&);

} // end namespace detail


//
// amida<T, C>
//

template<class T, std::size_t C>
class amida
{
public:
  typedef std::size_t                    size_type;
  typedef T                              value_type;
  typedef T&                             reference;
  typedef const T&                       const_reference;
  typedef detail::amida_pointer<T>       pointer;
  typedef const detail::amida_pointer<T> const_pointer;
  typedef amida_series_iterator<T>       series_iterator;
  typedef const amida_series_iterator<T> const_series_iterator;
  
  typedef detail::amida_node<T, boost::is_class<T>::value> node_type;

  // constructors
  explicit amida(size_type r = 0) :
    _array(0), _vacant(0), _series(0), _nodes(0), _links(0), _cuts(0), 
    _largest(0) { init(r); }
  amida(const amida<T>& a); // not implemented!
  ~amida() {}

  void init(size_type r);
  void clear() { init(series()); }

  size_type series() const { return _series; }
  size_type nodes() const { return _nodes; }
  size_type links() const { return _links; }
  size_type cuts() const { return _cuts; }
  size_type nodes_max() const { return _largest; }
  size_type capacity() const { return _array.capacity(); }
  double memory() const { return _array.memory(); }

  reference operator[](size_type n) { return _array[n]; }
  const_reference operator[](size_type n) const { return _array[n]; }

  pointer ptr(size_type n) { return pointer(&_array[n]); }
  const_pointer ptr(size_type n) const { return pointer(&_array[n]); }

  std::pair<series_iterator, series_iterator> series(size_type r);
  std::pair<const_series_iterator, const_series_iterator> series(size_type r) const;

  pointer
  insert_link(const T& t,
	      const series_iterator& prev0, const series_iterator& prev1,
	      const series_iterator& next0, const series_iterator& next1);
  pointer
  insert_link_next(const T& t,
		   const series_iterator& prev0, const series_iterator& prev1);
  pointer
  insert_link_prev(const T& t,
		   const series_iterator& next0, const series_iterator& next1);

  series_iterator
  insert_cut(const T& t,
	     const series_iterator& prev, const series_iterator& next);
  series_iterator
  insert_cut_next(const T& t, const series_iterator& prev);
  series_iterator
  insert_cut_prev(const T& t, const series_iterator& next);

  void erase(node_type* b);
  void erase(size_type n) { erase(&(_array[n])); }
  void erase(pointer p) { erase(static_cast<node_type*>(detail::addr(p))); }

  bool no_links(size_type r) const;

  inline friend void
  rotate_series(amida<T, C>& c, size_type first, size_type last, int direc) {
    // check parameter
    if (first > last || last > c.series())
      boost::throw_exception(std::invalid_argument("rotate_series() : invalid argument"));
    
    typedef amida<T, C>::node_type node_type;
    
    size_type s = c.series();
    size_type n = last - first;
    while (direc < 0) direc += n;
    std::vector<node_type> tmp(n);
    
    // rotate bottom
    for (size_type i = 0; i < n; ++i) tmp[i] = c._array[first + i];
    for (size_type i = 0; i < n; ++i)
      c._array[first + (direc + i) % n] = tmp[i];
    
    // rotate top
    for (size_type i = 0; i < n; ++i) tmp[i] = c._array[s + first + i];
    for (size_type i = 0; i < n; ++i)
      c._array[s + first + (direc + i) % n] = tmp[i];
    
    // adjust series[], next[], prev[]
    for (size_type i = 0; i < c._array.size(); ++i) {
      if (c._array[i].is_node()) {
	size_type t0 = c._array[i].series[0];
	if (t0 >= first && t0 < last)
	  c._array[i].series[0] = first + (direc - first + t0) % n;
	
	size_type t1 = c._array[i].series[1];
	if (t1 >= first && t1 < last)
	  c._array[i].series[1] = first + (direc - first + t1) % n;
	
	if (c._array[i].prev[0] != 0) {
	  size_type idx = index(c._array[i].prev[0], c);
	  if (idx >= first && idx < last)
	    c._array[i].prev[0] = &c._array[first + (direc - first + idx) % n];
	}
	if (c._array[i].prev[1] != 0) {
	  size_type idx = index(c._array[i].prev[1], c);
	  if (idx >= first && idx < last)
	    c._array[i].prev[1] = &c._array[first + (direc - first + idx) % n];
	}
	if (c._array[i].next[0] != 0) {
	  size_type idx = index(c._array[i].next[0], c);
	  if (idx >= s + first && idx < s + last)
	    c._array[i].next[0] =
	      &c._array[s + first + (direc - first + idx) % n];
	}
	if (c._array[i].next[1] != 0) {
	  size_type idx = index(c._array[i].next[1], c);
	  if (idx >= s + first && idx < s + last)
	    c._array[i].next[1] =
	      &c._array[s + first + (direc - first + idx) % n];
	}
      }
    }
  }

  inline friend size_type
  index(const void* p, const amida<T, C>& c) { 
    size_type r = std::numeric_limits<size_type>::max();
    for (int i = 0; i < c._array.size(); ++i) {
      if (&(c._array[i]) == p) {
	r = i;
	break; 
      }
    }
    return r;
  }
  inline friend size_type
  index(const pointer p, const amida<T, C>& c) {
    return index(detail::addr(p), c);
  }
  
  void output(std::ostream& os) const { output_serial(os); output_stack(os); }
  void output_serial(std::ostream& os = std::cout) const;
  void output_stack(std::ostream& os = std::cout) const;

  void save(alps::ODump& od) const;
  void load(alps::IDump& id);

protected:
  node_type* new_node();
  void delete_node(node_type* node);
  
private:
  std::deque<node_type> _array;
  node_type* _vacant; // pointer to the top of vacant node stack
  size_type _series;  // number of series
  size_type _nodes;   // number of nodes
  size_type _links;   // number of links
  size_type _cuts;    // number of cuts
  size_type _largest; // max number of nodes
};


namespace detail {

//
// struct amida_node_base
//

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


//
// amida_node<T, B>
//

// for non-intrinsic type

template<class T>
struct amida_node<T, true> : public T, public amida_node_base
{
  amida_node() : T(), amida_node_base() {}

  T& operator*() { return *this; }
  const T& operator*() const { return *this; }

  template<std::size_t C>
  void output(std::ostream& os, const amida<T, C>& c) const;

  template<std::size_t C>
  void save(alps::ODump& od, const amida<T, C>& c) const;
  template<std::size_t C>
  void load(alps::IDump& id, const amida<T, C>& c);
};

// for intrinsic types

template<class T>
struct amida_node<T, false> : amida_node_base
{
  amida_node() : amida_node_base(), _data() {}
  T _data;

  T& operator*() { return _data; }
  const T& operator*() const { return _data; }

  template<std::size_t C>
  void output(std::ostream& os, const amida<T, C>& c) const;

  template<std::size_t C>
  void save(alps::ODump& od, const amida<T, C>& c) const;
  template<std::size_t C>
  void load(alps::IDump& id, const amida<T, C>& c);
};


//
// amida_pointer<T>
//

template<class T>
class amida_pointer
{
public:
  typedef T value_type;
  typedef amida_node<T, boost::is_class<T>::value> node_type;

  friend amida_node_base* looper::detail::addr<T>(const amida_pointer<T>&);

  // constructors
  amida_pointer() : _node(0) {}
  amida_pointer(const amida_pointer<T>& s) : _node(s._node) {}
  amida_pointer(const node_type* n) : _node(const_cast<node_type*>(n)) {}

  T& operator*() { return _node->operator*(); }
  const T& operator*() const { return _node->operator*(); }
  node_type* operator->() { return _node; }
  const node_type* operator->() const { return _node; }

  bool operator==(const amida_pointer& i) const { return (_node == i._node); }
  bool operator!=(const amida_pointer& i) const { return (_node != i._node); }

protected:
  node_type* _node;
};

} // end namespace detail


//
// amida_series_iterator<T>
//

template<class T>
class amida_series_iterator : public detail::amida_pointer<T>
{
public:
  typedef std::bidirectional_iterator_tag iterator_category;
  typedef detail::amida_node<T, boost::is_class<T>::value> node_type;

  friend std::size_t looper::series<T>(const amida_series_iterator<T>& itr);
  friend std::size_t looper::leg<T>(const amida_series_iterator<T>& itr);
  friend void looper::jump<T>(amida_series_iterator<T>&);
  friend void looper::proceed<T>(amida_series_iterator<T>&, size_t);

  // constructors
  amida_series_iterator() : detail::amida_pointer<T>(), _leg(0), _series(0) {}
  amida_series_iterator(const amida_series_iterator<T>& s) :
    detail::amida_pointer<T>(s._node), _leg(s._leg), _series(s._series) {}
  amida_series_iterator(const detail::amida_pointer<T>& p, std::size_t l) :
    detail::amida_pointer<T>(p), _leg(l), _series(p->series[l]) {}
  amida_series_iterator(const node_type* n, std::size_t l) : 
    detail::amida_pointer<T>(n), _leg(l), _series(n->series[l]) {}

  void incr();
  void decr();
  amida_series_iterator& operator++();
  amida_series_iterator operator++(int);
  amida_series_iterator& operator--();
  amida_series_iterator operator--(int);
  
  T& operator*() { return detail::amida_pointer<T>::operator*(); }
  const T& operator*() const { return detail::amida_pointer<T>::operator*(); }
  node_type* operator->() { return detail::amida_pointer<T>::operator->(); }
  const node_type* operator->() const { return detail::amida_pointer<T>::operator->(); }

  bool operator==(const amida_series_iterator& i) const {
    return ((_node == i._node) && (_leg == i._leg)); }
  bool operator!=(const amida_series_iterator& i) const {
    return ((_node != i._node) || (_leg != i._leg)); }

  // fast access to the next and previous links
  detail::amida_pointer<T> next() {
    return static_cast<node_type*>(_node->_next[_leg]); }
  const detail::amida_pointer<T> next() const {
    return static_cast<node_type*>(_node->_next[_leg]); }

private:
  std::size_t _leg;
  std::size_t _series;
};


//
// member functions of amida<T, C>
//

template<class T, std::size_t C>
void amida<T, C>::init(std::size_t r) {
  // setup boundary nodes
  _series = r;
  _array.resize(2 * _series);
  for (std::size_t r = 0; r < _series; ++r) {
    std::size_t t = r + _series;
    _array[r].set_as_bottom(&_array[t], r);
    _array[t].set_as_top(&_array[r], r);
  }
  
  // setup stack of vacant nodes
  _largest = 2 * _series;
  _vacant = 0;

  // set current number of nodes, etc.
  _nodes = 2 * _series;
  _links = 0;
  _cuts = 0;
}

template<class T, std::size_t C>
inline std::pair<typename amida<T, C>::series_iterator, typename amida<T, C>::series_iterator>
amida<T, C>::series(std::size_t r) {
#ifdef LOOPER_DEBUG
  if (r >= _series)
    LOOPER_EXCEPTION("amida<>::series() : invalid range of r");
#endif
  return std::make_pair(series_iterator(&(_array[r]), 0),
			series_iterator(&(_array[r + _series]), 0));
}

template<class T, std::size_t C>
inline std::pair<typename amida<T, C>::const_series_iterator, typename amida<T, C>::const_series_iterator>
amida<T, C>::series(std::size_t r) const {
#ifdef LOOPER_DEBUG
  if (r >= _series)
    LOOPER_EXCEPTION("amida<>::series() : invalid range of r");
#endif
  return std::make_pair(series_iterator(&(_array[r]), 0),
			series_iterator(&(_array[r + _series]), 0));
}

template<class T, std::size_t C>
inline typename amida<T, C>::pointer
amida<T, C>::insert_link(const T& t,
  const series_iterator& prev0, const series_iterator& prev1,
  const series_iterator& next0, const series_iterator& next1) {
  // not optimized for intrinsic types

  // get a new node
  node_type* node = new_node();
  node->operator*() = t;
  
  node->series[0] = looper::series(prev0);
  node->series[1] = looper::series(prev1);
  
  // adjust link pointers
  node->prev[0] = detail::addr(prev0);
  node->prev[1] = detail::addr(prev1);
  detail::addr(prev0)->next[leg(prev0)] = node;
  detail::addr(prev1)->next[leg(prev1)] = node;
  
  node->next[0] = detail::addr(next0);
  node->next[1] = detail::addr(next1);
  detail::addr(next0)->prev[leg(next0)] = node;
  detail::addr(next1)->prev[leg(next1)] = node;
  
  ++_links;
  return pointer(node);
}

template<class T, std::size_t C>
inline typename amida<T, C>::pointer
amida<T, C>::insert_link_next(const T& t,
  const series_iterator& prev0, const series_iterator& prev1) {
  // not optimized for intrinsic types
  return insert_link(t, prev0, prev1,
		     ++series_iterator(prev0), ++series_iterator(prev1));
}

template<class T, std::size_t C>
inline typename amida<T, C>::pointer
amida<T, C>::insert_link_prev(const T& t,
  const series_iterator& next0, const series_iterator& next1) {
  // not optimized for intrinsic types
  return insert_link(t, --series_iterator(next0), --series_iterator(next1),
		     next0, next1);
}

template<class T, std::size_t C>
inline typename amida<T, C>::series_iterator
amida<T, C>::insert_cut(const T& t,
  const series_iterator& prev, const series_iterator& next) {
  // not optimized for intrinsic types

  // get a new node
  node_type* node = new_node();
  node->operator*() = t;
  
  node->series[0] = looper::series(prev);
  node->series[1] = looper::series(prev);
  
  // adjust link pointers
  node->prev[0] = detail::addr(prev);
  node->prev[1] = 0;
  detail::addr(prev)->next[leg(prev)] = node;
  
  node->next[0] = detail::addr(next);
  node->next[1] = detail::addr(next);
  detail::addr(next)->prev[leg(next)] = node;
  
  ++_cuts;
  return series_iterator(node, 0);
}

template<class T, std::size_t C>
inline typename amida<T, C>::series_iterator
amida<T, C>::insert_cut_next(const T& t, const series_iterator& prev) {
  // not optimized for intrinsic types
  return insert_cut(t, prev, ++series_iterator(prev));
}

template<class T, std::size_t C>
inline typename amida<T, C>::series_iterator
amida<T, C>::insert_cut_prev(const T& t, const series_iterator& next) {
  // not optimized for intrinsic types
  return insert_cut(t, --series_iterator(next), next);
}

template<class T, std::size_t C>
inline void
amida<T, C>::erase(node_type* node) {
  if (node->is_link()) {

    series_iterator prev0 = --series_iterator(node, 0);
    series_iterator prev1 = --series_iterator(node, 1);
    series_iterator next0 = ++series_iterator(node, 0);
    series_iterator next1 = ++series_iterator(node, 1);
    
    detail::addr(prev0)->next[leg(prev0)] = detail::addr(next0);
    detail::addr(prev1)->next[leg(prev1)] = detail::addr(next1);
    detail::addr(next0)->prev[leg(next0)] = detail::addr(prev0);
    detail::addr(next1)->prev[leg(next1)] = detail::addr(prev1);

    --_links;
  } else if (node->is_cut()) {

    series_iterator prev = --series_iterator(node, 0);
    series_iterator next = ++series_iterator(node, 0);
    
    detail::addr(prev)->next[leg(prev)] = detail::addr(next);
    detail::addr(next)->prev[leg(next)] = detail::addr(prev);

    --_cuts;
  } else {
    boost::throw_exception(std::logic_error("amida<>::erase() : couldn't erase node"));
  }
    
  delete_node(node);
}

template<class T, std::size_t C>
inline bool
amida<T, C>::no_links(std::size_t r) const {
#ifdef LOOPER_DEBUG
  if (r >= _series)
    LOOPER_EXCEPTION("amida<>::no_linkds() : invalid range of r");
#endif
  return (_array[r].next[0]->at_top());
}


template<class T, std::size_t C>
inline typename amida<T, C>::node_type*
amida<T, C>::new_node() {
  ++_nodes;
  if (_vacant == 0) {
    // stack expansion
    _array.resize(_array.size() + 1);
    node_type* tmp = &(_array[_array.size() - 1]);
    ++_largest;
    return tmp;
  } else {
    // pop from stack
    node_type* tmp = _vacant;
    _vacant = static_cast<node_type*>(_vacant->next[0]);
    return tmp;
  }
}

template<class T, std::size_t C>
inline void
amida<T, C>::delete_node(node_type* b) {
  // push to stack of vacant nodes
  b->set_as_vacant(_vacant);
  _vacant = b;
  --_nodes;
}

template<class T, std::size_t C>
void
amida<T, C>::save(alps::ODump& od) const {
  od << _largest << _series << _nodes << _links << _cuts
     << ((_vacant == 0) ? 0 : index(_vacant, *this) + 1);
  for (std::size_t i = 0; i < _largest; ++i) _array[i].save(od, *this);
}

template<class T, std::size_t C>
void
amida<T, C>::load(alps::IDump& id) {
  std::size_t vacant;
  id >> _largest >> _series >> _nodes >> _links >> _cuts >> vacant;
  _array.resize(_largest);
  _vacant  = static_cast<node_type*>(vacant == 0 ? 0 : &(_array[vacant - 1]));
  for (std::size_t i = 0; i < _largest; ++i) _array[i].load(id, *this);
}

template<class T, std::size_t C>
void
amida<T, C>::output_serial(std::ostream& os) const {
  os << "[amida Information]\n";
  for (std::size_t i = 0; i != _largest; ++i) {
    pointer p(this->ptr(i));
    series_iterator s0(p, 0);
    series_iterator s1(p, 1);
    if (!p->is_vacant()) {
      if (p->at_boundary()) {
	if (p->at_bottom()) {
  	  os << index(s0, *this) << '\t' << "root at " << looper::series(s0);
	  os << " :\t" << "next node is "
	     << index(++s0, *this);
	  os << "  ";
	  p->output(os, *this);
	  os << std::endl;
	} else {
	  os << index(s0, *this) << '\t' << "goal at " << looper::series(s0);
	  os << " :\t" << "prev node is " 
	     << index(--s0, *this);
	  os << "  ";
	  p->output(os, *this);
	  os << std::endl;
	}
      } else {
	os << index(s0, *this) << '\t' << "node connecting "
	   << looper::series(s0) << " and " << looper::series(s1);
	os << " :\t" << index((++s0)--, *this);
	os << ' ' << index((--s0)++, *this);
	os << ' ' << index((++s1)--, *this);
	os << ' ' << index((--s1)++, *this);
	os << '\t';
	p->output(os, *this);
	os << std::endl;
      }
    }
  }
}

template<class T, std::size_t C>
void
amida<T, C>::output_stack(std::ostream& os) const {
  os << "[amida Stack Information]\n";
  for (node_type* s = _vacant; s != 0; s = (node_type*)(s->next[0])) {
    os << index(s, *this) << "->";
  }
  os << "NULL\n";
}


//
// member functions of amida_node<T, B>
//

template<class T>
template<std::size_t C>
void
detail::amida_node<T, true>::output(std::ostream& os, const amida<T, C>& c) const {
  os << "node: " << index(this, c) << ", contents: ";
  T::output(os);
}

template<class T>
template<std::size_t C>
void
detail::amida_node<T, true>::save(alps::ODump& od, const amida<T, C>& c) const {
  od << series[0] << series[1]
     << ((next[0] == 0) ? 0 : index(next[0], c) + 1)
     << ((next[1] == 0) ? 0 : index(next[1], c) + 1)
     << ((prev[0] == 0) ? 0 : index(prev[0], c) + 1)
     << ((prev[1] == 0) ? 0 : index(prev[1], c) + 1);
  T::save(od);
}

template<class T>
template<std::size_t C>
void
detail::amida_node<T, true>::load(alps::IDump& id, const amida<T, C>& c) {
  std::size_t n0, n1, p0, p1;
  id >> series[0] >> series[1] >> n0 >> n1 >> p0 >> p1;
  next[0] = (n0 == 0 ? 0 : detail::addr(c.ptr(n0 - 1)));
  next[1] = (n1 == 0 ? 0 : detail::addr(c.ptr(n1 - 1)));
  prev[0] = (p0 == 0 ? 0 : detail::addr(c.ptr(p0 - 1)));
  prev[1] = (p1 == 0 ? 0 : detail::addr(c.ptr(p1 - 1)));
  T::load(id);
}

template<class T>
template<std::size_t C>
void
detail::amida_node<T, false>::output(std::ostream& os, const amida<T, C>& c) const {
  os << "node: " << index(this, c) << ", contents: " << _data;
}

template<class T>
template<std::size_t C>
void
detail::amida_node<T, false>::save(alps::ODump& od, const amida<T, C>& c) const {
  od << series[0] << series[1]
     << ((next[0] == 0) ? 0 : index(next[0], c) + 1)
     << ((next[1] == 0) ? 0 : index(next[1], c) + 1)
     << ((prev[0] == 0) ? 0 : index(prev[0], c) + 1)
     << ((prev[1] == 0) ? 0 : index(prev[1], c) + 1);
  od << _data;
}

template<class T>
template<std::size_t C>
void 
detail::amida_node<T, false>::load(alps::IDump& id, const amida<T, C>& c) {
  std::size_t n0, n1, p0, p1;
  id >> series[0] >> series[1] >> n0 >> n1 >> p0 >> p1;
  next[0] = (n0 == 0 ? 0 : detail::addr(c.ptr(n0 - 1)));
  next[1] = (n1 == 0 ? 0 : detail::addr(c.ptr(n1 - 1)));
  prev[0] = (p0 == 0 ? 0 : detail::addr(c.ptr(p0 - 1)));
  prev[1] = (p1 == 0 ? 0 : detail::addr(c.ptr(p1 - 1)));
  id >> _data;
}


//
// member functions of amida_series_iterator<T>
//

template<class T>
inline void
amida_series_iterator<T>::incr() {
#ifdef LOOPER_DEBUG
  if (_node == 0)
    LOOPER_EXCEPTION("amida_series_iterator::incr() : no node associated");
  if (_node->at_top())
    LOOPER_EXCEPTION("amida_series_iterator::incr() : at boundary");
#endif
  _node = static_cast<node_type*>(_node->next[_leg]);
  (_node->series[0] == _series) ? _leg = 0 : _leg = 1;
}

template<class T>
inline void
amida_series_iterator<T>::decr() {
#ifdef LOOPER_DEBUG
  if (_node == 0)
    LOOPER_EXCEPTION("amida_series_iterator::decr() : no node associated");
  if (_node->at_bottom())
    LOOPER_EXCEPTION("amida_series_iterator::decr() : at boundary");
#endif
  _node = static_cast<node_type*>(_node->prev[_leg]);
  (_node->series[0] == _series) ? _leg = 0 : _leg = 1;
}

template<class T>
inline amida_series_iterator<T>&
amida_series_iterator<T>::operator++() {
  incr(); return *this;
}

template<class T>
inline amida_series_iterator<T>
amida_series_iterator<T>::operator++(int) {
  amida_series_iterator tmp = *this;
  incr();
  return tmp;
}

template<class T>
inline amida_series_iterator<T>&
amida_series_iterator<T>::operator--() {
  decr(); return *this;
}

template<class T>
inline amida_series_iterator<T>
amida_series_iterator<T>::operator--(int) {
  amida_series_iterator tmp = *this;
  decr();
  return tmp;
}


//
// helper functions
//

template<class T>
inline detail::amida_node_base*
detail::addr(const detail::amida_pointer<T>& p) { return p._node; }

template<class T>
inline std::size_t 
series(const amida_series_iterator<T>& p) { return p._series; }
  
template<class T>
inline std::size_t
leg(const amida_series_iterator<T>& p) { return p._leg; }

template<class T>
inline void
jump(amida_series_iterator<T>& p) {
  p._leg ^= 1;
  p._series = p._node->series[p._leg];
}

template<class T>
inline void 
proceed(amida_series_iterator<T>& p, std::size_t d) {
  // d=0 for ++, 1 for --
  d ? p.decr() : p.incr(); 
} 


//
// I/O functions
//

template<class T, std::size_t C>
std::ostream& operator<<(std::ostream& os, const looper::amida<T, C>& c) {
  c.output_serial(os); c.output_stack(os); return os;
}

template<class T, std::size_t C>
alps::ODump& operator<<(alps::ODump& od, const looper::amida<T, C>& c) {
  c.save(od); return od;
}

template<class T, std::size_t C>
alps::IDump& operator>>(alps::IDump& id, looper::amida<T, C>& c) {
  c.load(id); return id;
}

} // end namespace looper

#endif // LOOPER_AMIDA_H
