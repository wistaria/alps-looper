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

#ifndef LOOPER_GRAPH_H_
#define LOOPER_GRAPH_H_

#include <alps/lattice.h>
#include <alps/model.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <utility>                 // std::pair, std::make_pair
#include <vector>                  // std::vector

namespace looper {

typedef alps::graph_name_t        graph_name_t;
typedef alps::dimension_t         dimension_t;
typedef boost::vertex_index_t     vertex_index_t;
typedef alps::vertex_type_t       vertex_type_t;
typedef alps::coordinate_t        coordinate_t;
typedef alps::parity_t            parity_t;
typedef boost::edge_index_t       edge_index_t;
typedef alps::edge_type_t         edge_type_t;
typedef alps::edge_vector_t       edge_vector_t;
typedef alps::boundary_crossing_t boundary_crossing_t;

// NOTE: We use adjacency_list with EdgeListS=vecS, which will be less
// efficient than that with EdgeListS=listS for edge insertion.
// However, the former allows us random access to edges, which is
// required for QMC in SSE representation.

typedef boost::adjacency_list<
  boost::vecS, boost::vecS, boost::undirectedS,
  boost::property<coordinate_t, alps::coordinate_type,
    boost::property<vertex_type_t, int > >,
  boost::property<edge_type_t, int,
    boost::property<edge_index_t, int,
      boost::property<boundary_crossing_t, alps::boundary_crossing,
        boost::property<edge_vector_t, alps::coordinate_type> > > >,
  boost::property<dimension_t, std::size_t,
    boost::property<graph_name_t, std::string > >,
  boost::vecS> graph_type;

typedef boost::adjacency_list<
  boost::vecS, boost::vecS, boost::undirectedS,
  boost::property<coordinate_t, alps::coordinate_type,
    boost::property<parity_t, int,
      boost::property<vertex_type_t, int> > >,
  boost::property<edge_type_t, int,
    boost::property<edge_index_t, int,
      boost::property<boundary_crossing_t, alps::boundary_crossing,
        boost::property<edge_vector_t, alps::coordinate_type> > > >,
  boost::property<dimension_t, std::size_t,
    boost::property<graph_name_t, std::string > >,
  boost::vecS> parity_graph_type;

template <class G>
struct graph_traits
{
  typedef G                                           graph_type;
  typedef std::pair<typename boost::graph_traits<graph_type>::vertex_iterator,
                    typename boost::graph_traits<graph_type>::vertex_iterator>
                                                      vertex_pair_type;
  typedef std::vector<std::vector<vertex_pair_type> > vp_type;
  typedef typename boost::property_traits<
    typename boost::property_map<graph_type, coordinate_t>::type>::value_type
                                                      coordinate_type;
};

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
    boost::graph_traits<unit_cell_type::graph_type>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(uc.graph());
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
  boost::get_property(g, graph_name_t()) = "simple hypercubic graph";
  alps::make_graph_from_lattice(g,
    lattice_type(cl, desc.extent().begin(), desc.extent().end(),
                 bc.begin(), bc.end()));
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline int
gauge(typename boost::graph_traits<
        boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
        vertex_descriptor vd,
      const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  BOOST_STATIC_ASSERT(alps::parity::white == 0);
  BOOST_STATIC_ASSERT(alps::parity::black == 1);
  return 1 - 2 * (int)boost::get(parity_t(), g, vd);
}

struct uniform
{
  template<class T0, class T1, class T2, class T3, class T4, class T5,
           class T6>
  static double
  value(typename boost::graph_traits<
          boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
          vertex_descriptor,
        const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>&) {
    return 1.;
  }
};

struct staggered
{
  template<class T0, class T1, class T2, class T3, class T4, class T5,
           class T6>
  static double
  value(typename boost::graph_traits<
          boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
          vertex_descriptor vd,
        const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g) {
    return 1. - 2. * boost::get(parity_t(), g, vd);
  }
};

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline unsigned int
site_index(typename boost::graph_traits<
             boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
             vertex_descriptor vd,
           const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  return boost::get(boost::vertex_index, g, vd);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline unsigned int
site_type(typename boost::graph_traits<
            boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
            vertex_descriptor vd,
          const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  return boost::get(alps::site_type_t(), g, vd);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline unsigned int
bond_index(typename boost::graph_traits<
             boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
             edge_descriptor ed,
           const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  return boost::get(boost::edge_index, g, ed);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline unsigned int
bond_type(typename boost::graph_traits<
            boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
            edge_descriptor ed,
          const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  return boost::get(alps::bond_type_t(), g, ed);
}


//
// class template virtual_mapping
//

// class template virtual_mapping
// for describing a mapping from a real site/bond to virtual ones

template<class G>
class virtual_mapping
{
public:
  typedef G                                           graph_type;
  typedef typename boost::graph_traits<graph_type>::vertex_descriptor
    vertex_descriptor;
  typedef typename boost::graph_traits<graph_type>::edge_descriptor
    edge_descriptor;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::edge_iterator
    edge_iterator;
  typedef std::pair<vertex_iterator, vertex_iterator> vertex_range_type;
  typedef std::pair<edge_iterator, edge_iterator>     edge_range_type;

  virtual_mapping() :
    vertex_map_(1, vertex_iterator()), edge_map_(1, edge_iterator()) {}

  vertex_range_type
  virtual_vertices(const graph_type& rg, const vertex_descriptor& rv) const {
    return std::make_pair(vertex_map_[boost::get(boost::vertex_index, rg, rv)],
      vertex_map_[boost::get(boost::vertex_index, rg, rv) + 1]);
  }

  edge_range_type
  virtual_edges(const graph_type& rg, const edge_descriptor& re) const {
    return std::make_pair(edge_map_[boost::get(boost::edge_index, rg, re)],
      edge_map_[boost::get(boost::edge_index, rg, re) + 1]);
  }

  void add_vertex(const graph_type& rg, const vertex_descriptor& rv,
                  const vertex_iterator& first, const vertex_iterator& last)
  {
    assert(boost::get(boost::vertex_index, rg, rv) == vertex_map_.size() - 1);
    if (vertex_map_.size() == 1) {
      vertex_map_.back() = first;
    } else {
      assert(first == vertex_map_.back());
    }
    vertex_map_.push_back(last);
  }

  void add_edge(const graph_type& rg, const edge_descriptor& re,
                const edge_iterator& first, const edge_iterator& last)
  {
    assert(boost::get(boost::edge_index, rg, re) == edge_map_.size() - 1);
    if (edge_map_.size() == 1) {
      edge_map_.back() = first;
    } else {
      assert(first == edge_map_.back());
    }
    edge_map_.push_back(last);
  }

  void clear()
  {
    vertex_map_.clear();
    edge_map_.clear();
    vertex_map_.push_back(vertex_iterator());
    edge_map_.push_back(edge_iterator());
  }

  bool operator==(const virtual_mapping& rhs) const
  { return vertex_map_ == rhs.vertex_map_ && edge_map_ == rhs.edge_map_; }
  bool operator!=(const virtual_mapping& rhs) { return !(*this == rhs); }

  void output(std::ostream& os, const graph_type& rg,
              const graph_type& vg) const
  {
    os << "[[vitual_mapping]]\n";
    os << "  number of vertex groups = " << boost::num_vertices(rg) << '\n';
    os << "  vertex mapping:\n";
    vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(rg); vi != vi_end; ++vi) {
      os << "    " << boost::get(boost::vertex_index, rg, *vi) << " -> ";
      vertex_range_type vr = virtual_vertices(rg, *vi);
      if (vr.first == vr.second) {
        os << "null\n";
      } else if (vr.first == boost::prior(vr.second)) {
        os << boost::get(boost::vertex_index, vg, *vr.first) << '\n';
      } else {
        os << '['
           << boost::get(boost::vertex_index, vg, *vr.first) << ','
           << boost::get(boost::vertex_index, vg, *boost::prior(vr.second))
           << "]\n";
      }
    }
    os << "  number of edge groups = " << boost::num_edges(rg) << '\n';
    os << "  edge mapping:\n";
    edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(rg); ei != ei_end; ++ei) {
      os << "    " << boost::get(boost::edge_index, rg, *ei) << " -> ";
      edge_range_type er = virtual_edges(rg, *ei);
      if (er.first == er.second) {
        os << "null\n";
      } else if (er.first == boost::prior(er.second)) {
        os << boost::get(boost::edge_index, vg, *er.first) << '\n';
      } else {
        os << '['
           << boost::get(boost::edge_index, vg, *er.first) << ','
           << boost::get(boost::edge_index, vg, *boost::prior(er.second))
           << "]\n";
      }
    }
  }

private:
  std::vector<vertex_iterator> vertex_map_;
  std::vector<edge_iterator> edge_map_;
};


//
// function generate_virtual_graph
//

template<class RealGraph, class M, class VirtualGraph>
inline void generate_virtual_graph(const RealGraph& rg,
                                   const M& model,
                                   VirtualGraph& vg,
                                   virtual_mapping<VirtualGraph>& vm)
{
  typedef RealGraph    rgraph_type;
  typedef VirtualGraph vgraph_type;

  vg.clear();
  vm.clear();

  // setup graph properties
  boost::get_property(vg, graph_name_t())
    = "virtual graph of " + boost::get_property(rg, graph_name_t());
  boost::get_property(vg, dimension_t())
    = boost::get_property(rg, dimension_t());

  // setup vertices
  {
    typename boost::graph_traits<rgraph_type>::vertex_iterator rvi, rvi_end;
    for (boost::tie(rvi, rvi_end) = boost::vertices(rg);
         rvi != rvi_end; ++rvi) {
      for (int i = 0; i < model.site(*rvi, rg).s().get_twice(); ++i) {
        // add vertices to virtual graph
        typename boost::graph_traits<vgraph_type>::vertex_descriptor
          vvd = boost::add_vertex(vg);
        
        // copy vertex properties
        alps::copy_property(vertex_type_t(), rg, *rvi, vg, vvd);
        alps::copy_property(coordinate_t(), rg, *rvi, vg, vvd);
        alps::copy_property(parity_t(), rg, *rvi, vg, vvd);
      }
    }
    // setup mapping
    typename boost::graph_traits<vgraph_type>::vertex_iterator
      vvi_first = boost::vertices(vg).first;
    typename boost::graph_traits<vgraph_type>::vertex_iterator
      vvi_last = vvi_first;
    for (boost::tie(rvi, rvi_end) = boost::vertices(rg);
         rvi != rvi_end; ++rvi) {
      vvi_last += model.site(*rvi, rg).s().get_twice();
      vm.add_vertex(rg, *rvi, vvi_first, vvi_last);
      vvi_first = vvi_last;
    }
  }

  // setup edges
  {
    typename boost::graph_traits<rgraph_type>::edge_iterator rei, rei_end;
    for (boost::tie(rei, rei_end) = boost::edges(rg);
	 rei != rei_end; ++rei) {
      typename boost::graph_traits<rgraph_type>::vertex_descriptor
        rs = boost::source(*rei, rg);
      typename boost::graph_traits<rgraph_type>::vertex_descriptor
        rt = boost::target(*rei, rg);
      typename boost::graph_traits<vgraph_type>::vertex_iterator
        vvsi, vvsi_end;
      for (boost::tie(vvsi, vvsi_end) = vm.virtual_vertices(rg, rs);
           vvsi != vvsi_end; ++vvsi) {
        typename boost::graph_traits<vgraph_type>::vertex_iterator
          vvti, vvti_end;
        for (boost::tie(vvti, vvti_end) = vm.virtual_vertices(rg, rt);
             vvti != vvti_end; ++vvti) {
          // add edges to virtual graph
          typename boost::graph_traits<vgraph_type>::edge_descriptor ved =
            boost::add_edge(*vvsi, *vvti, vg).first;
          
          // setup edge properties
          boost::put(boost::edge_index, vg, ved, boost::num_edges(vg) - 1);
          alps::copy_property(edge_type_t(), rg, *rei, vg, ved);
          alps::copy_property(edge_vector_t(), rg, *rei, vg, ved);
        }
      }
    }
    // setup mapping
    typename boost::graph_traits<vgraph_type>::edge_iterator
      vei_first = boost::edges(vg).first;
    typename boost::graph_traits<vgraph_type>::edge_iterator
      vei_last = vei_first;
    for (boost::tie(rei, rei_end) = boost::edges(rg);
	 rei != rei_end; ++rei) {
      typename boost::graph_traits<vgraph_type>::vertex_descriptor
	v0 = boost::source(*rei, rg);
      typename boost::graph_traits<vgraph_type>::vertex_descriptor
	v1 = boost::target(*rei, rg);
      vei_last += model.site(v0, rg).s().get_twice() *
	model.site(v1, rg).s().get_twice();
      vm.add_edge(rg, *rei, vei_first, vei_last);
      vei_first = vei_last;
    }
  }
}

namespace vg_detail {

template<class T>
struct spin_wrapper
{
  spin_wrapper(const T& v) : val_(v) {}
  const T& s() const { return val_; }
  const T& val_;
};

template<class T>
struct vector_spin_wrapper
{
  vector_spin_wrapper(const std::vector<T>& v) : vec_(v) {}
  template<class V, class G>
  spin_wrapper<T> site(const V& v, const G& g) const
  { return spin_wrapper<T>(vec_[boost::get(vertex_type_t(), g, v)]); }
  const std::vector<T>& vec_;
};

template<class T>
struct const_spin_wrapper
{
  const_spin_wrapper(const T& t) : t_(t) {}
  template<class V, class G>
  spin_wrapper<T> site(const V&, const G&) const
  { return spin_wrapper<T>(t_); }
  const T& t_;
};

}

template<class RealGraph, class IntType, class VirtualGraph>
inline void generate_virtual_graph(const RealGraph& rg,
                                   const alps::half_integer<IntType>& s,
                                   VirtualGraph& vg,
                                   virtual_mapping<VirtualGraph>& vm)
{
  generate_virtual_graph(rg,
    vg_detail::const_spin_wrapper<alps::half_integer<IntType> >(s),
    vg, vm);
}

template<class RealGraph, class IntType, class VirtualGraph>
inline void generate_virtual_graph(const RealGraph& rg,
  const std::vector<alps::half_integer<IntType> >& v,
  VirtualGraph& vg, virtual_mapping<VirtualGraph>& vm)
{
  generate_virtual_graph(rg,
    vg_detail::vector_spin_wrapper<alps::half_integer<IntType> >(v),
    vg, vm);
}

} // end namespace looper

#endif // LOOPER_GRAPH_H
