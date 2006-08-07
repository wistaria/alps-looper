/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_LATTICE_H_
#define LOOPER_LATTICE_H_

#include "integer_range.h"
#include <alps/lattice.h>
#include <alps/math.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <utility>                 // std::pair, std::make_pair
#include <vector>                  // std::vector

namespace looper {

using alps::vertex_index_t;
using alps::site_index_t;
using alps::vertex_type_t;
using alps::site_type_t;
using alps::coordinate_t;
using alps::parity_t;

using alps::edge_index_t;
using alps::bond_index_t;
using alps::edge_type_t;
using alps::bond_type_t;
using alps::edge_vector_t;
using alps::bond_vector_t;
using alps::edge_vector_relative_t;
using alps::bond_vector_relative_t;

using alps::graph_name_t;
using alps::dimension_t;

using alps::graph_traits;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
  // vertex property
  boost::property<vertex_type_t, alps::type_type,
  boost::property<coordinate_t, alps::coordinate_type,
  boost::property<parity_t, int> > >,
  // edge property
  boost::property<edge_index_t, unsigned int,
  boost::property<edge_type_t, alps::type_type,
  boost::property<edge_vector_t, alps::coordinate_type,
  boost::property<edge_vector_relative_t, alps::coordinate_type> > > >,
  // graph property
  boost::property<dimension_t, std::size_t,
  boost::property<graph_name_t, std::string > >,
  boost::vecS> default_graph_type;

struct real_vertex_index_t {typedef boost::vertex_property_tag kind; };
typedef real_vertex_index_t real_site_index_t;
struct real_edge_index_t {typedef boost::edge_property_tag kind; };
typedef real_edge_index_t real_bond_index_t;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
  // vertex property
  boost::property<real_vertex_index_t, int>,
  // edge property
  boost::property<edge_index_t, unsigned int,
  boost::property<real_edge_index_t, int> >,
  // graph property
  boost::no_property,
  boost::vecS> virtual_graph_type;

template<class D = std::size_t, class S = D, class E = std::vector<S> >
class hypercubic_graph_generator
{
public:
  typedef D dimension_type;
  typedef S size_type;
  typedef E extent_type;

  hypercubic_graph_generator() : ext_() {}
  explicit hypercubic_graph_generator(dimension_type dim, size_type size)
    : ext_(dim, size) {}
  template<class Container>
  hypercubic_graph_generator(const Container& ext) : ext_(ext.size()) {
    std::copy(ext.begin(), ext.end(), ext_.begin());
  }
  dimension_type dimension() const { return ext_.size(); }
  size_type length(dimension_type d) const { return ext_[d]; }
  const extent_type& extent() const { return ext_; }

private:
  extent_type ext_;
};


template<class D, class S, class E, class G>
void generate_graph(G& g, const hypercubic_graph_generator<D, S, E>& desc)
{
  std::size_t dim(desc.dimension());

  typedef alps::hypercubic_lattice<
    alps::coordinate_lattice<alps::simple_lattice<alps::GraphUnitCell> > >
    lattice_type;
  typedef alps::lattice_traits<lattice_type>::unit_cell_type unit_cell_type;

  // unit cell
  unit_cell_type uc(dim);
  uc.add_vertex(0, unit_cell_type::coordinate_type(dim, 0.0));
  alps::lattice_traits<lattice_type>::offset_type so(dim, 0), to(dim, 0);
  for (std::size_t d = 0; d < dim; ++d) {
    to[d] = 1; uc.add_edge(d, 1, so, 1, to); to[d] = 0;
  }

  // edge_vectors
  if (alps::has_property<edge_vector_t,unit_cell_type::graph_type>
        ::edge_property) {
    std::vector<double> b(dim, 0.0);
    int d = 0;
    graph_traits<unit_cell_type::graph_type>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(uc.graph());
         ei != ei_end; ++ei, ++d) {
      b[d] = 1.; alps::get_or_default(edge_vector_t(),uc.graph(),
                   alps::coordinate_type())[*ei] = b; b[d] = 0.;
    }
  }

  // basis vectors
  lattice_type::parent_lattice_type cl(uc);
  std::vector<double> b(dim, 0.0);
  for (std::size_t d = 0; d < dim; ++d) {
    b[d] = 1.0; cl.add_basis_vector(b); b[d] = 0.0;
  }

  // boundary condition
  std::vector<std::string> bc(dim, "periodic");
  for (std::size_t d = 0; d < dim; ++d)
    if (desc.extent()[d] <= 2) bc[d] = "open";

  // generate graph
  g.clear();
  get_property(g, graph_name_t()) = "simple hypercubic graph";
  alps::make_graph_from_lattice(g,
    lattice_type(cl, desc.extent().begin(), desc.extent().end(),
                 bc.begin(), bc.end()));
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline int
gauge(const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g,
      typename graph_traits<
        boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
          vertex_descriptor vd)
{ return 2 * get(parity_t(), g, vd) - 1; }


//
// class template virtual_mapping
//
// for describing a mapping from a real vertex/edge to virtual ones

template<typename RG, typename VG>
class virtual_mapping
{
public:
  typedef RG real_graph_type;
  typedef typename graph_traits<real_graph_type>::vertex_descriptor
    real_vertex_descriptor;
  typedef typename graph_traits<real_graph_type>::edge_descriptor
    real_edge_descriptor;
  typedef typename graph_traits<real_graph_type>::vertex_iterator
    real_vertex_iterator;
  typedef typename graph_traits<real_graph_type>::edge_iterator
    real_edge_iterator;

  typedef VG virtual_graph_type;
  typedef typename graph_traits<virtual_graph_type>::vertex_iterator
    virtual_vertex_iterator;
  typedef typename graph_traits<virtual_graph_type>::edge_iterator
    virtual_edge_iterator;

  typedef std::pair<virtual_vertex_iterator, virtual_vertex_iterator>
    virtual_vertex_range_type;
  typedef std::pair<virtual_edge_iterator, virtual_edge_iterator>
    virtual_edge_range_type;

  virtual_mapping() :
    vertex_map_(1, virtual_vertex_iterator()),
    edge_map_(1, virtual_edge_iterator()),
    v2e_map_(1, virtual_edge_iterator()),
    v2e_offset_(0), max_vv_(0) {}

  virtual_vertex_range_type
  virtual_vertices(const real_graph_type& rg,
                   const real_vertex_descriptor& rv) const {
    return std::make_pair(vertex_map_[get(vertex_index_t(), rg, rv)],
                          vertex_map_[get(vertex_index_t(), rg, rv) + 1]);
  }

  virtual_edge_range_type
  virtual_edges(const real_graph_type& rg,
                const real_edge_descriptor& re) const {
    return std::make_pair(edge_map_[get(edge_index_t(), rg, re)],
                          edge_map_[get(edge_index_t(), rg, re) + 1]);
  }

  virtual_edge_range_type
  virtual_edges(const real_graph_type& rg,
                const real_vertex_descriptor& rv) const {
    return std::make_pair(v2e_map_[get(vertex_index_t(), rg, rv)],
                          v2e_map_[get(vertex_index_t(), rg, rv) + 1]);
  }

  void add_vertices(const real_graph_type& rg,
                    const real_vertex_descriptor& rv,
                    const virtual_vertex_iterator& first,
                    const virtual_vertex_iterator& last)
  {
    if (get(vertex_index_t(), rg, rv) != vertex_map_.size() - 1)
      boost::throw_exception(std::invalid_argument(
        "virtual_mapping<G>::add_vertices()"));
    if (vertex_map_.size() == 1) {
      vertex_map_.back() = first;
    } else {
      assert(first == vertex_map_.back());
    }
    vertex_map_.push_back(last);
    max_vv_ = std::max(max_vv_, static_cast<int>(last - first));
  }

  void add_edges(const real_graph_type& rg,
                 const real_edge_descriptor& re,
                 const virtual_edge_iterator& first,
                 const virtual_edge_iterator& last)
  {
    if (get(edge_index_t(), rg, re) != edge_map_.size() - 1)
      boost::throw_exception(std::invalid_argument(
        "virtual_mapping<G>::add_edges()"));
    if (edge_map_.size() == 1) {
      edge_map_.back() = first;
    } else {
      assert(first == edge_map_.back());
    }
    edge_map_.push_back(last);
  }

  void add_v2edges(const real_graph_type& rg,
                   const real_vertex_descriptor& rv,
                   const virtual_edge_iterator& first,
                   const virtual_edge_iterator& last)
  {
    if (get(vertex_index_t(), rg, rv) != v2e_map_.size() - 1)
      boost::throw_exception(std::invalid_argument(
        "virtual_mapping<G>::add_v2edges()"));
    if (v2e_map_.size() == 1) {
      v2e_map_.back() = first;
    } else {
      assert(first == v2e_map_.back());
    }
    v2e_map_.push_back(last);
  }

  void clear()
  {
    vertex_map_.clear();
    vertex_map_.push_back(virtual_vertex_iterator());
    edge_map_.clear();
    edge_map_.push_back(virtual_edge_iterator());
    v2e_map_.clear();
    v2e_map_.push_back(virtual_edge_iterator());
    v2e_offset_ = 0;
    max_vv_ = 0;
  }

  void set_v2edge_type_offset(int t) { v2e_offset_ = t; }
  int v2edge_type_offset() const { return v2e_offset_; }

  int max_virtual_vertices() const { return max_vv_; }

  bool operator==(const virtual_mapping& rhs) const
  { return vertex_map_ == rhs.vertex_map_ && edge_map_ == rhs.edge_map_ &&
      v2e_map_ == rhs.v2e_map_; }
  bool operator!=(const virtual_mapping& rhs) { return !(*this == rhs); }

  void output(std::ostream& os, const real_graph_type& rg,
              const virtual_graph_type& vg) const;

private:
  std::vector<virtual_vertex_iterator> vertex_map_;
  std::vector<virtual_edge_iterator> edge_map_;
  std::vector<virtual_edge_iterator> v2e_map_;
  int v2e_offset_;
  int max_vv_;
};


//
// class template virtual_lattice
//

template<typename RG, typename VG = virtual_graph_type>
class virtual_lattice {
public:
  typedef RG real_graph_type;
  typedef typename graph_traits<real_graph_type>::vertex_descriptor
    real_vertex_descriptor;
  typedef typename graph_traits<real_graph_type>::edge_descriptor
    real_edge_descriptor;
  typedef typename graph_traits<real_graph_type>::vertex_iterator
    real_vertex_iterator;
  typedef typename graph_traits<real_graph_type>::edge_iterator
    real_edge_iterator;
  typedef real_vertex_descriptor real_site_descriptor;
  typedef real_edge_descriptor   real_bond_descriptor;
  typedef real_vertex_iterator   real_site_iterator;
  typedef real_edge_iterator     real_bond_iterator;

  typedef VG virtual_graph_type;
  typedef typename graph_traits<virtual_graph_type>::vertex_iterator
    virtual_vertex_iterator;
  typedef typename graph_traits<virtual_graph_type>::edge_iterator
    virtual_edge_iterator;
  typedef typename graph_traits<virtual_graph_type>::vertex_descriptor
    virtual_vertex_descriptor;
  typedef typename graph_traits<virtual_graph_type>::edge_descriptor
    virtual_edge_descriptor;
  typedef virtual_vertex_iterator   virtual_site_iterator;
  typedef virtual_edge_iterator     virtual_bond_iterator;
  typedef virtual_vertex_descriptor virtual_site_descriptor;
  typedef virtual_edge_descriptor    virtual_bond_descriptor;

  typedef virtual_mapping<real_graph_type, virtual_graph_type> mapping_type;
  typedef typename mapping_type::virtual_vertex_range_type
    virtual_vertex_range_type;
  typedef typename mapping_type::virtual_edge_range_type
    virtual_edge_range_type;
  typedef virtual_vertex_range_type virtual_site_range_type;
  typedef virtual_edge_range_type   virtual_bond_range_type;

  virtual_lattice(const real_graph_type& rg) : rgraph_(rg) {}
  template<class M>
  virtual_lattice(const real_graph_type& rg, const M& model,
                  bool has_d_term = false) :
    rgraph_(rg)
  { reinitialize(model, has_d_term); }

  // real vertex -> virtual vertex
  virtual_vertex_range_type
  virtual_vertices(const real_vertex_descriptor& rv) const
  { return mapping_.virtual_vertices(rgraph_, rv); }

  // real edge -> virtual edge
  virtual_edge_range_type
  virtual_edges(const real_edge_descriptor& re) const
  { return mapping_.virtual_edges(rgraph_, re); }

  // real vertex -> virtual edge (for high spins with D term)
  virtual_edge_range_type
  virtual_edges(const real_vertex_descriptor& rv) const
  { return mapping_.virtual_edges(rgraph_, rv); }

  void clear() { vgraph_.clear(); mapping_.clear(); }

  template<class M>
  void reinitialize(const M& model, bool has_d_term = false);

  const real_graph_type& rgraph() const { return rgraph_; }
  const virtual_graph_type& vgraph() const { return vgraph_; }
  const mapping_type& mapping() const { return mapping_; }

  void print_mapping(std::ostream& os) const {
    mapping_.output(os, rgraph_, vgraph_);
  }

private:
  real_graph_type const& rgraph_;
  virtual_graph_type vgraph_;
  mapping_type mapping_;
  integer_range<int> vtype_range_, etype_range_;
};

template<typename RG, typename VG>
typename graph_traits<VG>::vertices_size_type
num_vvertices(const virtual_lattice<RG, VG>& vl)
{ return num_vertices(vl.vgraph()); }

template<typename RG, typename VG>
typename graph_traits<VG>::edges_size_type
num_vedges(const virtual_lattice<RG, VG>& vl)
{ return num_edges(vl.vgraph()); }

template<typename RG, typename VG>
std::pair<typename graph_traits<VG>::vertex_iterator,
          typename graph_traits<VG>::vertex_iterator>
vvertices(const virtual_lattice<RG, VG>& vl) { return vertices(vl.vgraph()); }

template<typename RG, typename VG>
std::pair<typename graph_traits<VG>::edge_iterator,
          typename graph_traits<VG>::edge_iterator>
vedges(const virtual_lattice<RG, VG>& vl) { return edges(vl.vgraph()); }

template<typename RG, typename VG>
typename virtual_lattice<RG, VG>::virtual_vertex_range_type
virtual_vertices(const virtual_lattice<RG, VG>& vl,
                 const typename graph_traits<RG>::vertex_descriptor& rv)
{ return vl.virtual_vertices(rv); }

template<typename RG, typename VG>
typename virtual_lattice<RG, VG>::virtual_edge_range_type
virtual_edges(const virtual_lattice<RG, VG>& vl,
              const typename graph_traits<RG>::edge_descriptor& re)
{ return vl.virtual_edges(re); }

template<typename RG, typename VG>
typename virtual_lattice<RG, VG>::virtual_edge_range_type
virtual_edges(const virtual_lattice<RG, VG>& vl,
              const typename graph_traits<RG>::vertex_descriptor& rv)
{ return vl.virtual_edges(rv); }


template<typename RG, typename VG>
typename graph_traits<VG>::sites_size_type
num_vsites(const virtual_lattice<RG, VG>& vl)
{ return num_sites(vl.vgraph()); }

template<typename RG, typename VG>
typename graph_traits<VG>::bonds_size_type
num_vbonds(const virtual_lattice<RG, VG>& vl)
{ return num_bonds(vl.vgraph()); }

template<typename RG, typename VG>
std::pair<typename graph_traits<VG>::site_iterator,
          typename graph_traits<VG>::site_iterator>
vsites(const virtual_lattice<RG, VG>& vl) { return sites(vl.vgraph()); }

template<typename RG, typename VG>
std::pair<typename graph_traits<VG>::bond_iterator,
          typename graph_traits<VG>::bond_iterator>
vbonds(const virtual_lattice<RG, VG>& vl) { return bonds(vl.vgraph()); }

template<typename RG, typename VG>
typename virtual_lattice<RG, VG>::virtual_site_range_type
virtual_sites(const virtual_lattice<RG, VG>& vl,
              const typename graph_traits<RG>::site_descriptor& rv)
{ return vl.virtual_vertices(rv); }

template<typename RG, typename VG>
typename virtual_lattice<RG, VG>::virtual_bond_range_type
virtual_bonds(const virtual_lattice<RG, VG>& vl,
              const typename graph_traits<RG>::bond_descriptor& re)
{ return vl.virtual_edges(re); }

template<typename RG, typename VG>
typename virtual_lattice<RG, VG>::virtual_bond_range_type
virtual_bonds(const virtual_lattice<RG, VG>& vl,
              const typename graph_traits<RG>::site_descriptor& rv)
{ return vl.virtual_edges(rv); }


template<typename RG, typename VG>
typename graph_traits<VG>::site_descriptor
vsource(typename graph_traits<VG>::bond_descriptor b,
        const virtual_lattice<RG, VG>& vl)
{ return source(b, vl.vgraph()); }

template<typename RG, typename VG>
typename graph_traits<VG>::site_descriptor
vsource(typename graph_traits<VG>::bonds_size_type i,
        const virtual_lattice<RG, VG>& vl)
{ return source(bond(i, vl.vgraph()), vl.vgraph()); }

template<typename RG, typename VG>
typename graph_traits<VG>::site_descriptor
vtarget(typename graph_traits<VG>::bond_descriptor b,
        const virtual_lattice<RG, VG>& vl)
{ return target(b, vl.vgraph()); }

template<typename RG, typename VG>
typename graph_traits<VG>::site_descriptor
vtarget(typename graph_traits<VG>::bonds_size_type i,
        const virtual_lattice<RG, VG>& vl)
{ return target(bond(i, vl.vgraph()), vl.vgraph()); }

template<typename RG, typename VG>
inline bool
is_real_edge(const virtual_lattice<RG, VG>& vl,
             typename virtual_lattice<RG, VG>::virtual_edge_descriptor ved)
{ return (get(real_edge_index_t(), vl.vgraph(), ved) >= 0); }

template<typename RG, typename VG>
inline bool
is_real_bond(const virtual_lattice<RG, VG>& vl,
             typename virtual_lattice<RG, VG>::virtual_edge_descriptor ved)
{ return is_real_edge(vl, ved); }

template<typename RG, typename VG>
inline typename graph_traits<RG>::vertex_descriptor
rvertex(virtual_lattice<RG, VG> const& vl, unsigned int s)
{
  return *(vertices(vl.rgraph()).first +
           get(real_vertex_index_t(), vl.vgraph(), s));
}

template<typename RG, typename VG>
inline typename graph_traits<RG>::site_descriptor
rsite(virtual_lattice<RG, VG> const& vl, unsigned int s)
{ return rvertex(vl, s); }

template<typename RG, typename VG>
inline typename graph_traits<RG>::edge_descriptor
redge(virtual_lattice<RG, VG> const& vl,
      typename virtual_lattice<RG, VG>::virtual_edge_descriptor ved)
{
  return *(edges(vl.rgraph()).first +
           get(real_edge_index_t(), vl.vgraph(), ved));
}

template<typename RG, typename VG>
inline typename graph_traits<RG>::edge_descriptor
redge(virtual_lattice<RG, VG> const& vl, unsigned int b)
{ return redge(vl, *(edges(vl.vgraph()).first + b)); }

template<typename RG, typename VG>
inline typename graph_traits<RG>::bond_descriptor
rbond(virtual_lattice<RG, VG> const& vl,
      typename virtual_lattice<RG, VG>::virtual_bond_descriptor ved)
{ return redge(vl, ved); }

template<typename RG, typename VG>
inline typename graph_traits<RG>::bond_descriptor
rbond(virtual_lattice<RG, VG> const& vl, unsigned int b)
{ return redge(vl, b); }


template<typename RG, typename VG>
inline int
gauge(virtual_lattice<RG, VG> const& vl,
      typename virtual_lattice<RG, VG>::virtual_vertex_descriptor vvd)
{ return 2 * get(parity_t(), vl.rgraph(), rvertex(vl, vvd)) - 1; }


template<typename RG, typename VG>
int max_virtual_vertices(const virtual_lattice<RG, VG>& vl)
{ return vl.mapping().max_virtual_vertices(); }

template<typename RG, typename VG>
int max_virtual_sites(const virtual_lattice<RG, VG>& vl)
{ return vl.mapping().max_virtual_vertices(); }

template<typename RG, typename VG>
std::ostream& operator<<(std::ostream& os, const virtual_lattice<RG, VG>& vl)
{
  os << vl.vgraph();
  return os;
}

} // end namespace looper

namespace looper {

///////////////////////////////////////////////////////////////////////////////
//
// Implementations
//

//
// class template virtual_mapping
//

template<typename RG, typename VG>
void virtual_mapping<RG, VG>::output(std::ostream& os, RG const& rg,
                                     VG const& vg) const
{
  os << "[[vitual_mapping]]\n";
  os << "  number of vertex groups = " << num_vertices(rg) << '\n';
  os << "  vertex mapping:\n";
  virtual_vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(rg); vi != vi_end; ++vi) {
    os << "    " << get(vertex_index_t(), rg, *vi) << " -> ";
    virtual_vertex_range_type vr = virtual_vertices(rg, *vi);
    if (vr.first == vr.second) {
      os << "null\n";
    } else if (vr.first == boost::prior(vr.second)) {
      os << get(vertex_index_t(), vg, *vr.first) << '\n';
    } else {
      os << '['
         << get(vertex_index_t(), vg, *vr.first) << ','
         << get(vertex_index_t(), vg, *boost::prior(vr.second))
         << "]\n";
    }
  }
  os << "  number of edge groups = " << num_edges(rg) << '\n';
  os << "  edge mapping:\n";
  real_edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(rg); ei != ei_end; ++ei) {
    os << "    " << get(edge_index_t(), rg, *ei) << " -> ";
    virtual_edge_range_type er = virtual_edges(rg, *ei);
    if (er.first == er.second) {
      os << "null\n";
    } else if (er.first == boost::prior(er.second)) {
      os << get(edge_index_t(), vg, *er.first) << '\n';
    } else {
      os << '['
         << get(edge_index_t(), vg, *er.first) << ','
         << get(edge_index_t(), vg, *boost::prior(er.second))
         << "]\n";
    }
  }
  os << "  vertex2edge mapping:\n";
  for (boost::tie(vi, vi_end) = vertices(rg); vi != vi_end; ++vi) {
    os << "    " << get(vertex_index_t(), rg, *vi) << " -> ";
    virtual_edge_range_type er = virtual_edges(rg, *vi);
    if (er.first == er.second) {
      os << "null\n";
    } else if (er.first == boost::prior(er.second)) {
      os << get(edge_index_t(), vg, *er.first) << '\n';
    } else {
      os << '['
         << get(edge_index_t(), vg, *er.first) << ','
         << get(edge_index_t(), vg, *boost::prior(er.second))
         << "]\n";
    }
  }
}


//
// class template virtual_lattice
//

template<typename RG, typename VG>
template<class M>
void
virtual_lattice<RG, VG>::reinitialize(const M& model, bool has_d_term)
{
  vgraph_.clear();
  mapping_.clear();

  // setup v2edge_type_offset
  real_vertex_iterator rvi, rvi_end;
  int tmin = 0;
  for (boost::tie(rvi, rvi_end) = vertices(rgraph_); rvi != rvi_end; ++rvi)
    tmin = std::min(tmin, int(get(vertex_type_t(), rgraph_, *rvi)));
  real_edge_iterator rei, rei_end;
  int tmax = 0;
  for (boost::tie(rei, rei_end) = edges(rgraph_); rei != rei_end; ++rei)
    tmax = std::max(tmax, int(get(edge_type_t(), rgraph_, *rei)));
  mapping_.set_v2edge_type_offset(tmax-tmin+1);

  // add vertices to virtual graph
  for (boost::tie(rvi, rvi_end) = vertices(rgraph_); rvi != rvi_end; ++rvi)
    for (int i = 0; i < model.site(*rvi, rgraph_).s.get_twice(); ++i) {
      virtual_vertex_descriptor vvd = add_vertex(vgraph_);
      put(real_vertex_index_t(), vgraph_, vvd,
          get(vertex_index_t(), rgraph_, *rvi));
    }

  // setup vertex mapping
  virtual_vertex_iterator vvi_first = vertices(vgraph_).first;
  virtual_vertex_iterator vvi_last = vvi_first;
  for (boost::tie(rvi, rvi_end) = vertices(rgraph_); rvi != rvi_end; ++rvi) {
    vvi_last += model.site(*rvi, rgraph_).s.get_twice();
    mapping_.add_vertices(rgraph_, *rvi, vvi_first, vvi_last);
    vvi_first = vvi_last;
  }

  // add edges to virtual graph
  for (boost::tie(rei, rei_end) = edges(rgraph_); rei != rei_end; ++rei) {
    real_vertex_descriptor rs = source(*rei, rgraph_);
    real_vertex_descriptor rt = target(*rei, rgraph_);
    virtual_vertex_iterator vvsi, vvsi_end;
    for (boost::tie(vvsi, vvsi_end) = mapping_.virtual_vertices(rgraph_, rs);
         vvsi != vvsi_end; ++vvsi) {
      virtual_vertex_iterator vvti, vvti_end;
      for (boost::tie(vvti, vvti_end) = mapping_.virtual_vertices(rgraph_, rt);
           vvti != vvti_end; ++vvti) {
        virtual_edge_descriptor ved = add_edge(*vvsi, *vvti, vgraph_).first;
        put(edge_index_t(), vgraph_, ved, num_edges(vgraph_) - 1);
        put(real_edge_index_t(), vgraph_, ved,
            get(edge_index_t(), rgraph_, *rei));
      }
    }
  }

  // add `in-real-vertex' edges to virtual graph
  if (has_d_term)
    for (boost::tie(rvi, rvi_end) = vertices(rgraph_); rvi != rvi_end; ++rvi) {
      virtual_vertex_iterator vvsi, vvsi_end;
      for (boost::tie(vvsi, vvsi_end) =
             mapping_.virtual_vertices(rgraph_, *rvi); vvsi != vvsi_end;
           ++vvsi)
        for (virtual_vertex_iterator vvti = boost::next(vvsi);
             vvti != vvsi_end; ++vvti) {
          virtual_edge_descriptor ved = add_edge(*vvsi, *vvti, vgraph_).first;
          put(edge_index_t(), vgraph_, ved, num_edges(vgraph_) - 1);
          put(real_edge_index_t(), vgraph_, ved,
              get(vertex_index_t(), rgraph_, *rvi) - 1);
        }
    }

  // setup edge and v2edge mapping
  virtual_edge_iterator vei_first = edges(vgraph_).first;
  virtual_edge_iterator vei_last = vei_first;
  for (boost::tie(rei, rei_end) = edges(rgraph_); rei != rei_end; ++rei) {
    vei_last += model.site(source(*rei, rgraph_), rgraph_).s.get_twice() *
      model.site(target(*rei, rgraph_), rgraph_).s.get_twice();
    mapping_.add_edges(rgraph_, *rei, vei_first, vei_last);
    vei_first = vei_last;
  }
  for (boost::tie(rvi, rvi_end) = vertices(rgraph_); rvi != rvi_end; ++rvi) {
    if (has_d_term)
      vei_last += model.site(*rvi, rgraph_).s.get_twice() *
        (model.site(*rvi, rgraph_).s.get_twice() - 1) / 2;
    mapping_.add_v2edges(rgraph_, *rvi, vei_first, vei_last);
    vei_first = vei_last;
  }

  // set range of types
  vtype_range_ = 0;
  for (boost::tie(rvi, rvi_end) = vertices(rgraph_); rvi != rvi_end; ++rvi)
    vtype_range_.include(get(vertex_index_t(), rgraph_, *rvi));
  etype_range_ = 0;
  for (boost::tie(rei, rei_end) = edges(rgraph_); rei != rei_end; ++rei)
    etype_range_.include(get(edge_index_t(), rgraph_, *rei));
}

} // end namespace looper

#endif // LOOPER_LATTICE_H
