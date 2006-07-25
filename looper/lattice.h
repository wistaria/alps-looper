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
using alps::boundary_crossing_t;

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
  boost::property<boundary_crossing_t, alps::boundary_crossing,
  boost::property<edge_vector_t, alps::coordinate_type,
  boost::property<edge_vector_relative_t, alps::coordinate_type> > > > >,
  // graph property
  boost::property<dimension_t, std::size_t,
  boost::property<graph_name_t, std::string > >,
  boost::vecS> graph_type;


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

template<class G>
class virtual_mapping
{
public:
  typedef G graph_type;
  typedef typename graph_traits<graph_type>::vertex_descriptor
    vertex_descriptor;
  typedef typename graph_traits<graph_type>::edge_descriptor
    edge_descriptor;
  typedef typename graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename graph_traits<graph_type>::edge_iterator
    edge_iterator;
  typedef std::pair<vertex_iterator, vertex_iterator>
    vertex_range_type;
  typedef std::pair<edge_iterator, edge_iterator>
    edge_range_type;

  virtual_mapping()
    : vertex_map_(1, vertex_iterator()), edge_map_(1, edge_iterator()),
      v2e_map_(1, edge_iterator()), v2e_offset_(0), max_vv_(0) {}

  vertex_range_type
  virtual_vertices(int s) const {
    return std::make_pair(vertex_map_[s], vertex_map_[s+1]);
  }
  vertex_range_type
  virtual_vertices(const graph_type& rg, const vertex_descriptor& rv) const {
    return virtual_vertices(get(vertex_index_t(), rg, rv));
  }

  edge_range_type
  virtual_edges(const graph_type& rg, const edge_descriptor& re) const {
    return std::make_pair(edge_map_[get(edge_index_t(), rg, re)],
      edge_map_[get(edge_index_t(), rg, re) + 1]);
  }

  edge_range_type
  virtual_edges(const graph_type& rg, const vertex_descriptor& rv) const {
    return std::make_pair(v2e_map_[get(vertex_index_t(), rg, rv)],
      v2e_map_[get(vertex_index_t(), rg, rv) + 1]);
  }

  void add_vertices(const graph_type& rg, const vertex_descriptor& rv,
                    const vertex_iterator& first,
                    const vertex_iterator& last)
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

  void add_edges(const graph_type& rg, const edge_descriptor& re,
                 const edge_iterator& first, const edge_iterator& last)
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

  void add_v2edges(const graph_type& rg, const vertex_descriptor& rv,
                   const edge_iterator& first, const edge_iterator& last)
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
    vertex_map_.push_back(vertex_iterator());
    edge_map_.clear();
    edge_map_.push_back(edge_iterator());
    v2e_map_.clear();
    v2e_map_.push_back(edge_iterator());
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

  void output(std::ostream& os, const graph_type& rg,
              const graph_type& vg) const;

private:
  std::vector<vertex_iterator> vertex_map_;
  std::vector<edge_iterator> edge_map_;
  std::vector<edge_iterator> v2e_map_;
  int v2e_offset_;
  int max_vv_;
};


//
// class template virtual_lattice
//

template<class GRAPH>
class virtual_lattice {
public:
  typedef GRAPH                       graph_type;
  typedef virtual_mapping<graph_type> mapping_type;

  typedef typename graph_traits<graph_type>::vertex_descriptor
    vertex_descriptor;
  typedef typename graph_traits<graph_type>::edge_descriptor
    edge_descriptor;
  typedef typename graph_traits<graph_type>::adjacency_iterator
     adjacency_iterator;
  typedef typename graph_traits<graph_type>::out_edge_iterator
     out_edge_iterator;
  typedef typename graph_traits<graph_type>::in_edge_iterator
     in_edge_iterator;
  typedef typename graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename graph_traits<graph_type>::edge_iterator
    edge_iterator;

  typedef typename graph_traits<graph_type>::directed_category
    directed_category;
  typedef typename graph_traits<graph_type>::edge_parallel_category
    edge_parallel_category;
  typedef typename graph_traits<graph_type>::traversal_category
    traversal_category;

  typedef typename graph_traits<graph_type>::vertices_size_type
    vertices_size_type;
  typedef typename graph_traits<graph_type>::edges_size_type
    edges_size_type;
  typedef typename graph_traits<graph_type>::degree_size_type
    degree_size_type;

  virtual_lattice() {}
  template<class M>
  virtual_lattice(const graph_type& rg, const M& model,
                  bool has_d_term = false)
  { generate(rg, model, has_d_term); }

  // real vertex -> virtual vertex
  std::pair<vertex_iterator, vertex_iterator>
  virtual_vertices(int s) const
  { return mapping_.virtual_vertices(s); }
  std::pair<vertex_iterator, vertex_iterator>
  virtual_vertices(const graph_type& rg, const vertex_descriptor& rv) const
  { return mapping_.virtual_vertices(rg, rv); }

  // real edge -> virtual edge
  std::pair<edge_iterator, edge_iterator>
  virtual_edges(const graph_type& rg, const edge_descriptor& re) const
  { return mapping_.virtual_edges(rg, re); }

  // real vertex -> virtual edge (high spin with D term)
  std::pair<edge_iterator, edge_iterator>
  virtual_edges(const graph_type& rg, const vertex_descriptor& rv) const
  { return mapping_.virtual_edges(rg, rv); }

  void clear() { graph_.clear(); mapping_.clear(); }

  template<class M>
  void generate(const graph_type& rg, const M& model, bool has_d_term = false);

  const graph_type& graph() const { return graph_; }
  const mapping_type& mapping() const { return mapping_; }

  void print_mapping(std::ostream& os, const graph_type& rg) const {
    mapping_.output(os, rg, graph_);
  }

private:
  graph_type graph_;
  mapping_type mapping_;
};

template<class G>
typename graph_traits<virtual_lattice<G> >::vertices_size_type
num_vertices(const virtual_lattice<G>& vl)
{ return num_vertices(vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::edges_size_type
num_edges(const virtual_lattice<G>& vl)
{ return num_edges(vl.graph()); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::vertex_iterator,
          typename graph_traits<virtual_lattice<G> >::vertex_iterator>
vertices(const virtual_lattice<G>& vl)
{ return vertices(vl.graph()); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::edge_iterator,
          typename graph_traits<virtual_lattice<G> >::edge_iterator>
edges(const virtual_lattice<G>& vl)
{ return edges(vl.graph()); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::vertex_iterator,
          typename graph_traits<virtual_lattice<G> >::vertex_iterator>
virtual_vertices(const virtual_lattice<G>& vl, const G& rg,
                 const typename graph_traits<G>::vertex_descriptor& rv)
{ return vl.virtual_vertices(rg, rv); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::edge_iterator,
          typename graph_traits<virtual_lattice<G> >::edge_iterator>
virtual_edges(const virtual_lattice<G>& vl, const G& rg,
              const typename graph_traits<G>::edge_descriptor& re)
{ return vl.virtual_edges(rg, re); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::edge_iterator,
          typename graph_traits<virtual_lattice<G> >::edge_iterator>
virtual_edges(const virtual_lattice<G>& vl, const G& rg,
              const typename graph_traits<G>::vertex_descriptor& rv)
{ return vl.virtual_edges(rg, rv); }


template<class G>
typename graph_traits<virtual_lattice<G> >::sites_size_type
num_sites(const virtual_lattice<G>& vl)
{ return num_sites(vl.graph()); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::site_iterator,
          typename graph_traits<virtual_lattice<G> >::site_iterator>
sites(const virtual_lattice<G>& vl)
{ return sites(vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::site_descriptor
site(typename alps::graph_traits<virtual_lattice<G> >::sites_size_type i,
     const virtual_lattice<G>& vl)
{ return site(i, vl.graph()); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::bond_iterator,
          typename graph_traits<virtual_lattice<G> >::bond_iterator>
bonds(const virtual_lattice<G>& vl)
{ return bonds(vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::bond_descriptor
bond(typename alps::graph_traits<virtual_lattice<G> >::bonds_size_type i,
     const virtual_lattice<G>& vl)
{ return bond(i, vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::bonds_size_type
num_bonds(const virtual_lattice<G>& vl)
{ return num_edges(vl.graph()); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::site_iterator,
          typename graph_traits<virtual_lattice<G> >::site_iterator>
virtual_sites(const virtual_lattice<G>& vl, int s)
{ return vl.virtual_vertices(s); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::site_iterator,
          typename graph_traits<virtual_lattice<G> >::site_iterator>
virtual_sites(const virtual_lattice<G>& vl, const G& rg,
              const typename graph_traits<G>::site_descriptor& rv)
{ return vl.virtual_vertices(rg, rv); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::bond_iterator,
          typename graph_traits<virtual_lattice<G> >::bond_iterator>
virtual_bonds(const virtual_lattice<G>& vl, const G& rg,
              const typename graph_traits<G>::bond_descriptor& re)
{ return vl.virtual_edges(rg, re); }

template<class G>
std::pair<typename graph_traits<virtual_lattice<G> >::bond_iterator,
          typename graph_traits<virtual_lattice<G> >::bond_iterator>
virtual_bonds(const virtual_lattice<G>& vl, const G& rg,
              const typename graph_traits<G>::site_descriptor& rv)
{ return vl.virtual_edges(rg, rv); }

template<class G>
typename graph_traits<virtual_lattice<G> >::site_descriptor
vsource(typename alps::graph_traits<virtual_lattice<G> >::bond_descriptor b,
        const virtual_lattice<G>& vl)
{ return source(b, vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::site_descriptor
vsource(typename alps::graph_traits<virtual_lattice<G> >::bonds_size_type i,
        const virtual_lattice<G>& vl)
{ return source(bond(i, vl.graph()), vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::site_descriptor
vtarget(typename alps::graph_traits<virtual_lattice<G> >::bond_descriptor b,
        const virtual_lattice<G>& vl)
{ return target(b, vl.graph()); }

template<class G>
typename graph_traits<virtual_lattice<G> >::site_descriptor
vtarget(typename alps::graph_traits<virtual_lattice<G> >::bonds_size_type i,
        const virtual_lattice<G>& vl)
{ return target(bond(i, vl.graph()), vl.graph()); }

template<class G>
int gauge(const virtual_lattice<G>& vl,
          const typename virtual_lattice<G>::vertex_descriptor& vd)
{ return gauge(vl.graph(), vd); }

template<class G>
int max_virtual_vertices(const virtual_lattice<G>& vl)
{ return vl.mapping().max_virtual_vertices(); }

template<class G>
int max_virtual_sites(const virtual_lattice<G>& vl)
{ return vl.mapping().max_virtual_vertices(); }

template<class G>
std::ostream& operator<<(std::ostream& os, const virtual_lattice<G>& vl)
{
  os << vl.graph();
  return os;
}

//
// virtual_lattice_adaptor
//

template<class G>
class virtual_lattice_adaptor {
public:
  typedef G graph_type;
  typedef virtual_lattice<graph_type> virtual_lattice_type;

  virtual_lattice_adaptor(graph_type& rg) : rgraph_(rg) {}
  template<class MP>
  virtual_lattice_adaptor(graph_type& rg, const MP& mp)
    : rgraph_(rg), vlat_(rg, mp, mp.has_d_term()) {}

  template<class MP>
  void init(const MP& mp) { vlat_.generate(rgraph_, mp, mp.has_d_term()); }

  const graph_type& rgraph() const { return rgraph_; }
  const graph_type& vgraph() const { return vlat_.graph(); }
  const virtual_lattice_type& vlattice() const { return vlat_; }

private:
  graph_type& rgraph_;
  virtual_lattice_type vlat_;
};

} // end namespace looper

namespace looper {

////////////////////////////////////////////////////////////////////////////////
//
// Implementations
//

//
// class template virtual_mapping
//

template<class G>
void virtual_mapping<G>::output(std::ostream& os, const graph_type& rg,
  const graph_type& vg) const
{
  os << "[[vitual_mapping]]\n";
  os << "  number of vertex groups = " << num_vertices(rg) << '\n';
  os << "  vertex mapping:\n";
  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(rg); vi != vi_end; ++vi) {
    os << "    " << get(vertex_index_t(), rg, *vi) << " -> ";
    vertex_range_type vr = virtual_vertices(rg, *vi);
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
  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(rg); ei != ei_end; ++ei) {
    os << "    " << get(edge_index_t(), rg, *ei) << " -> ";
    edge_range_type er = virtual_edges(rg, *ei);
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
    edge_range_type er = virtual_edges(rg, *vi);
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

template<class G>
template<class M>
inline void
virtual_lattice<G>::generate(const G& rg, const M& model, bool has_d_term)
{
  graph_.clear();
  mapping_.clear();

  // setup graph properties
  get_property(graph_, graph_name_t()) = "virtual graph of " +
    get_property(rg, graph_name_t());
  get_property(graph_, dimension_t()) = get_property(rg, dimension_t());

  // setup v2edge_type_offset
  vertex_iterator rvi, rvi_end;
  int tmin = 0;
  for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi)
    tmin = std::min(tmin, int(get(vertex_type_t(), rg, *rvi)));
  edge_iterator rei, rei_end;
  int tmax = 0;
  for (boost::tie(rei, rei_end) = edges(rg); rei != rei_end; ++rei)
    tmax = std::max(tmax, int(get(edge_type_t(), rg, *rei)));
  mapping_.set_v2edge_type_offset(tmax-tmin+1);

  // add vertices to virtual graph
  for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi) {
    for (int i = 0; i < model.site(*rvi, rg).s.get_twice(); ++i) {
      vertex_descriptor vvd = add_vertex(graph_);
      copy_property(vertex_type_t(), rg, *rvi, graph_, vvd);
      copy_property(coordinate_t(), rg, *rvi, graph_, vvd);
      copy_property(parity_t(), rg, *rvi, graph_, vvd);
    }
  }

  // setup vertex mapping
  vertex_iterator vvi_first = vertices(graph_).first;
  vertex_iterator vvi_last = vvi_first;
  for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi) {
    vvi_last += model.site(*rvi, rg).s.get_twice();
    mapping_.add_vertices(rg, *rvi, vvi_first, vvi_last);
    vvi_first = vvi_last;
  }

  // add edges to virtual graph
  for (boost::tie(rei, rei_end) = edges(rg); rei != rei_end; ++rei) {
    vertex_descriptor rs = source(*rei, rg);
    vertex_descriptor rt = target(*rei, rg);
    vertex_iterator vvsi, vvsi_end;
    for (boost::tie(vvsi, vvsi_end) = mapping_.virtual_vertices(rg, rs);
         vvsi != vvsi_end; ++vvsi) {
      vertex_iterator vvti, vvti_end;
      for (boost::tie(vvti, vvti_end) = mapping_.virtual_vertices(rg, rt);
           vvti != vvti_end; ++vvti) {
        edge_descriptor ved = add_edge(*vvsi, *vvti, graph_).first;
        put(edge_index_t(), graph_, ved, num_edges(graph_) - 1);
        copy_property(edge_type_t(), rg, *rei, graph_, ved);
        copy_property(boundary_crossing_t(), rg, *rei, graph_, ved);
        copy_property(edge_vector_t(), rg, *rei, graph_, ved);
        copy_property(edge_vector_relative_t(), rg, *rei, graph_, ved);
      }
    }
  }

  // add `in-real-vertex' edges to virtual graph
  if (has_d_term) {
    int dim = alps::get_or_default(dimension_t(), rg, int(0));
    alps::coordinate_type vec(dim, 0);
    for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi) {
      int t = mapping_.v2edge_type_offset() +
        get(vertex_type_t(), rg, *rvi);
      vertex_iterator vvsi, vvsi_end;
      for (boost::tie(vvsi, vvsi_end) = mapping_.virtual_vertices(rg, *rvi);
           vvsi != vvsi_end; ++vvsi) {
        for (vertex_iterator vvti = boost::next(vvsi); vvti != vvsi_end;
             ++vvti) {
          edge_descriptor ved = add_edge(*vvsi, *vvti, graph_).first;
          put(edge_index_t(), graph_, ved, num_edges(graph_) - 1);
          put(edge_type_t(), graph_, ved, t);
          put(boundary_crossing_t(), graph_, ved,
                     alps::boundary_crossing());
          put(edge_vector_t(), graph_, ved, vec);
          put(edge_vector_relative_t(), graph_, ved, vec);
        }
      }
    }
  }

  // setup edge and v2edge mapping
  edge_iterator vei_first = edges(graph_).first;
  edge_iterator vei_last = vei_first;
  for (boost::tie(rei, rei_end) = edges(rg); rei != rei_end; ++rei) {
    vei_last += model.site(source(*rei, rg), rg).s.get_twice() *
      model.site(target(*rei, rg), rg).s.get_twice();
    mapping_.add_edges(rg, *rei, vei_first, vei_last);
    vei_first = vei_last;
  }
  for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi) {
    if (has_d_term) {
      vei_last += model.site(*rvi, rg).s.get_twice() *
        (model.site(*rvi, rg).s.get_twice() - 1) / 2;
    }
    mapping_.add_v2edges(rg, *rvi, vei_first, vei_last);
    vei_first = vei_last;
  }
}

} // end namespace looper

#endif // LOOPER_LATTICE_H
