// $Id: worldline.h 398 2003-10-09 10:33:05Z wistaria $

// looper-3 : C++ library for continuous-time loop algorithm
//
// Copyright (C) 2001,2002  Synge Todo <wistaria@itp.phys.ethz.ch>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

#ifndef LOOPER_WORLDLINE_H
#define LOOPER_WORLDLINE_H

#include <looper/amida.h>
#include <looper/node.h>
#include <looper/unionfind.h>
#include <alps/osiris.h>

namespace looper {

//
// class template WorldLine
// 

template<bool HasCTime = false>
class WorldLine
{
public:
  typedef detail::Node<HasCTime>                     node_type;
  typedef typename Amida<node_type>::pointer         node_pointer;
  typedef typename Amida<node_type>::series_iterator series_iterator;

  static const bool has_ctime = HasCTime;

  // constructors & destructor
  WorldLine() : _config() {}
  template<class GraphT>
  WorldLine(const GraphT& vg, int c = 0) : _config() { init(vg, c); }

  // initialize
  template<class GraphT>
  void init(const GraphT& vg, int c = 0)
  {
    assert(c == 0 || c == 1);
    int n = boost::num_vertices(vg);
    _config.init(n);
    for (std::size_t s = 0; s < n; ++s) {
      bottom(s)->set_time(0);
      bottom(s)->set_conf(c);
      top(s)->set_time(1);
      top(s)->set_conf(c);
    }
  }
  template<class GraphT, class RNG>
  void init(const GraphT& vg, RNG& rng, double p = 0.5)
  {
    init(vg);
    for (std::size_t s = 0; s < sites(); ++s) {
      if (rng() < p) {
	bottom(s)->set_conf(0);
	top(s)->set_conf(0);
      } else {
	bottom(s)->set_conf(1);
	top(s)->set_conf(1);
      }
    }
  }
  
  // size inquiry
  std::size_t sites() const { return _config.series(); }
  std::size_t links() const { return _config.links(); }
  std::size_t cuts() const { return _config.cuts(); }
  std::size_t nodes() const { return _config.nodes(); }
  std::size_t nodes_max() const { return _config.nodes_max(); }
  
  double memory() const { return _config.memory(); }

  // generate iterators, pointers, etc.
  series_iterator bottom(std::size_t r) { return _config.series(r).first; }
  const series_iterator bottom(std::size_t r) const {
    return _config.series(r).first;
  }
  series_iterator top(std::size_t r) { return _config.series(r).second; }
  const series_iterator top(std::size_t r) const {
    return _config.series(r).second;
  }
  node_pointer node(std::size_t i) { return _config.ptr(i); }
  const node_pointer node(std::size_t i) const { return _config.ptr(i); }

  std::pair<series_iterator, series_iterator>
  insert(const series_iterator& curr0, const series_iterator& curr1,
	 const series_iterator& next0, const series_iterator& next1,
	 double t)
  {
    // insert to list
    node_type k;
    k.set_time(t);
    node_pointer itr_new = _config.insert_link(k, curr0, curr1, next0, next1);
    return std::make_pair(series_iterator(itr_new, 0),
			  series_iterator(itr_new, 1));
  }

  void erase(node_pointer link) { _config.erase(link); }
  
  void save(alps::ODump& od) const { _config.save(od); }
  void load(alps::IDump& id) { _config.load(id); }
  
private:
  Amida<node_type> _config;
};

} // end namespace looper


#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<bool HasCTime>
alps::ODump& operator<<(alps::ODump& od, const WorldLine<HasCTime>& wline) {
  wline.save(od);
  return od;
}

template<bool HasCTime>
alps::IDump& operator>>(alps::IDump& id, WorldLine<HasCTime>& wline) {
  wline.load(id);
  return id;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_WORLDLINE_H
