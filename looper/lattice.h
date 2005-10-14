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

#ifndef LOOPER_LATTICE_H_
#define LOOPER_LATTICE_H_

#include <alps/lattice.h>
#include <alps/math.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <utility>                 // std::pair, std::make_pair
#include <vector>                  // std::vector

namespace looper {

using alps::site_index_t;
using alps::site_type_t;
using alps::coordinate_t;
struct parity_t { typedef boost::vertex_property_tag kind; };

using alps::bond_index_t;
using alps::bond_type_t;
using alps::bond_vector_t;
using alps::bond_vector_relative_t;
using alps::boundary_crossing_t;

using alps::graph_name_t;
using alps::dimension_t;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
  // vertex property
  boost::property<site_type_t, alps::type_type,
  boost::property<coordinate_t, alps::coordinate_type,
  boost::property<parity_t, int> > >,
  // edge property
  boost::property<bond_index_t, unsigned int,
  boost::property<bond_type_t, alps::type_type,
  boost::property<boundary_crossing_t, alps::boundary_crossing,
  boost::property<bond_vector_t, alps::coordinate_type,
  boost::property<bond_vector_relative_t, alps::coordinate_type> > > > >,
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

  // bond_vectors
  if (alps::has_property<bond_vector_t,unit_cell_type::graph_type>
        ::bond_property) {
    std::vector<double> b(dim, 0.0);
    int d = 0;
    alps::graph_traits<unit_cell_type::graph_type>::bond_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = bonds(uc.graph());
         ei != ei_end; ++ei, ++d) {
      b[d] = 1.; alps::get_or_default(bond_vector_t(),uc.graph(),
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
gauge(const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g,
      typename alps::graph_traits<
        boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
          site_descriptor vd)
{
  return boost::get(parity_t(), g, vd);
}

// struct uniform
// {
//   template<class T0, class T1, class T2, class T3, class T4, class T5,
//            class T6>
//   static double
//   value(typename alps::graph_traits<
//           boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
//           site_descriptor,
//         const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>&) {
//     return 1.;
//   }
// };

// struct staggered
// {
//   template<class T0, class T1, class T2, class T3, class T4, class T5,
//            class T6>
//   static int
//   value(typename alps::graph_traits<
//           boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
//           site_descriptor vd,
//         const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g) {
//     return boost::get(parity_t(), g, vd);
//   }
// };

// template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
// inline unsigned int
// site_index(typename alps::graph_traits<
//              boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
//              site_descriptor vd,
//            const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
// {
//   return boost::get(site_index_t(), g, vd);
// }

// template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
// inline unsigned int
// site_type(typename alps::graph_traits<
//             boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
//             site_descriptor vd,
//           const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
// {
//   return boost::get(alps::site_type_t(), g, vd);
// }

// template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
// inline unsigned int
// bond_index(typename alps::graph_traits<
//              boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
//              bond_descriptor ed,
//            const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
// {
//   return boost::get(boost::bond_index, g, ed);
// }

// template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
// inline unsigned int
// bond_type(typename alps::graph_traits<
//             boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
//             bond_descriptor ed,
//           const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
// {
//   return boost::get(alps::bond_type_t(), g, ed);
// }


//
// class template virtual_mapping
//
// for describing a mapping from a real site/bond to virtual ones
template<class G>
class virtual_mapping
{
public:
  typedef G graph_type;
  typedef typename alps::graph_traits<graph_type>::site_descriptor
    site_descriptor;
  typedef typename alps::graph_traits<graph_type>::bond_descriptor
    bond_descriptor;
  typedef typename alps::graph_traits<graph_type>::site_iterator
    site_iterator;
  typedef typename alps::graph_traits<graph_type>::bond_iterator
    bond_iterator;
  typedef std::pair<site_iterator, site_iterator>
    site_range_type;
  typedef std::pair<bond_iterator, bond_iterator>
    bond_range_type;

  virtual_mapping()
    : site_map_(1, site_iterator()), bond_map_(1, bond_iterator()),
      v2e_map_(1, bond_iterator())
  {}

  site_range_type
  virtual_sites(const graph_type& rg, const site_descriptor& rv) const {
    return std::make_pair(site_map_[boost::get(site_index_t(), rg, rv)],
      site_map_[boost::get(site_index_t(), rg, rv) + 1]);
  }

  bond_range_type
  virtual_bonds(const graph_type& rg, const bond_descriptor& re) const {
    return std::make_pair(bond_map_[boost::get(bond_index_t(), rg, re)],
      bond_map_[boost::get(bond_index_t(), rg, re) + 1]);
  }

  bond_range_type
  virtual_bonds(const graph_type& rg, const site_descriptor& rv) const {
    return std::make_pair(v2e_map_[boost::get(site_index_t(), rg, rv)],
      v2e_map_[boost::get(site_index_t(), rg, rv) + 1]);
  }

  void add_sites(const graph_type& rg, const site_descriptor& rv,
                    const site_iterator& first, const site_iterator& last)
  {
    assert(boost::get(site_index_t(), rg, rv) == site_map_.size() - 1);
    if (site_map_.size() == 1) {
      site_map_.back() = first;
    } else {
      assert(first == site_map_.back());
    }
    site_map_.push_back(last);
  }

  void add_bonds(const graph_type& rg, const bond_descriptor& re,
                 const bond_iterator& first, const bond_iterator& last)
  {
    assert(boost::get(bond_index_t(), rg, re) == bond_map_.size() - 1);
    if (bond_map_.size() == 1) {
      bond_map_.back() = first;
    } else {
      assert(first == bond_map_.back());
    }
    bond_map_.push_back(last);
  }

  void add_v2bonds(const graph_type& rg, const site_descriptor& rv,
                   const bond_iterator& first, const bond_iterator& last)
  {
    assert(boost::get(site_index_t(), rg, rv) == v2e_map_.size() - 1);
    if (v2e_map_.size() == 1) {
      v2e_map_.back() = first;
    } else {
      assert(first == v2e_map_.back());
    }
    v2e_map_.push_back(last);
  }

  void clear()
  {
    site_map_.clear();
    site_map_.push_back(site_iterator());
    bond_map_.clear();
    bond_map_.push_back(bond_iterator());
    v2e_map_.clear();
    v2e_map_.push_back(bond_iterator());
    v2e_offset_ = 0;
  }

  void set_v2bond_type_offset(int t) { v2e_offset_ = t; }
  int v2bond_type_offset() const { return v2e_offset_; }

  bool operator==(const virtual_mapping& rhs) const
  { return site_map_ == rhs.site_map_ && bond_map_ == rhs.bond_map_ &&
      v2e_map_ == rhs.v2e_map_; }
  bool operator!=(const virtual_mapping& rhs) { return !(*this == rhs); }

  void output(std::ostream& os, const graph_type& rg,
              const graph_type& vg) const
  {
    os << "[[vitual_mapping]]\n";
    os << "  number of site groups = " << num_sites(rg) << '\n';
    os << "  site mapping:\n";
    site_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = sites(rg); vi != vi_end; ++vi) {
      os << "    " << boost::get(site_index_t(), rg, *vi) << " -> ";
      site_range_type vr = virtual_sites(rg, *vi);
      if (vr.first == vr.second) {
        os << "null\n";
      } else if (vr.first == boost::prior(vr.second)) {
        os << boost::get(site_index_t(), vg, *vr.first) << '\n';
      } else {
        os << '['
           << boost::get(site_index_t(), vg, *vr.first) << ','
           << boost::get(site_index_t(), vg, *boost::prior(vr.second))
           << "]\n";
      }
    }
    os << "  number of bond groups = " << num_bonds(rg) << '\n';
    os << "  bond mapping:\n";
    bond_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = alps::bonds(rg); ei != ei_end; ++ei) {
      os << "    " << boost::get(bond_index_t(), rg, *ei) << " -> ";
      bond_range_type er = virtual_bonds(rg, *ei);
      if (er.first == er.second) {
        os << "null\n";
      } else if (er.first == boost::prior(er.second)) {
        os << boost::get(bond_index_t(), vg, *er.first) << '\n';
      } else {
        os << '['
           << boost::get(bond_index_t(), vg, *er.first) << ','
           << boost::get(bond_index_t(), vg, *boost::prior(er.second))
           << "]\n";
      }
    }
    os << "  site2bond mapping:\n";
    for (boost::tie(vi, vi_end) = sites(rg); vi != vi_end; ++vi) {
      os << "    " << boost::get(site_index_t(), rg, *vi) << " -> ";
      bond_range_type er = virtual_bonds(rg, *vi);
      if (er.first == er.second) {
        os << "null\n";
      } else if (er.first == boost::prior(er.second)) {
        os << boost::get(bond_index_t(), vg, *er.first) << '\n';
      } else {
        os << '['
           << boost::get(bond_index_t(), vg, *er.first) << ','
           << boost::get(bond_index_t(), vg, *boost::prior(er.second))
           << "]\n";
      }
    }
  }

private:
  std::vector<site_iterator> site_map_;
  std::vector<bond_iterator> bond_map_;
  std::vector<bond_iterator> v2e_map_;
  int v2e_offset_;
};


//
// function generate_virtual_lattice
//

template<class G, class M, class W>
inline void generate_virtual_lattice(G& vg, virtual_mapping<G>& vm,
  const G& rg, const M& model, const W& weight)
{
  // using boost::get_property;
  // using alps::copy_property;

  typedef G graph_type;
  typedef typename alps::graph_traits<G>::site_iterator site_iterator;
  typedef typename alps::graph_traits<G>::site_descriptor site_descriptor;
  typedef typename alps::graph_traits<G>::bond_iterator bond_iterator;
  typedef typename alps::graph_traits<G>::bond_descriptor bond_descriptor;

  vg.clear();
  vm.clear();

  // setup graph properties
  get_property(vg, graph_name_t()) = "virtual graph of " +
    get_property(rg, graph_name_t());
  get_property(vg, dimension_t()) = get_property(rg, dimension_t());

  // setup v2bond_type_offset
  site_iterator rvi, rvi_end;
  int tmin = 0;
  for (boost::tie(rvi, rvi_end) = sites(rg); rvi != rvi_end; ++rvi)
    tmin = std::min(tmin, int(boost::get(site_type_t(), rg, *rvi)));
  bond_iterator rei, rei_end;
  int tmax = 0;
  for (boost::tie(rei, rei_end) = bonds(rg); rei != rei_end; ++rei)
    tmax = std::max(tmax, int(boost::get(bond_type_t(), rg, *rei)));
  vm.set_v2bond_type_offset(tmax-tmin+1);

  // add sites to virtual graph
  for (boost::tie(rvi, rvi_end) = sites(rg); rvi != rvi_end; ++rvi) {
    for (int i = 0; i < model.site(*rvi, rg).s().get_twice(); ++i) {
      site_descriptor vvd = add_vertex(vg);
      copy_property(site_type_t(), rg, *rvi, vg, vvd);
      copy_property(coordinate_t(), rg, *rvi, vg, vvd);
      copy_property(parity_t(), rg, *rvi, vg, vvd);
    }
  }

  // setup site mapping
  site_iterator vvi_first = sites(vg).first;
  site_iterator vvi_last = vvi_first;
  for (boost::tie(rvi, rvi_end) = sites(rg); rvi != rvi_end; ++rvi) {
    vvi_last += model.site(*rvi, rg).s().get_twice();
    vm.add_sites(rg, *rvi, vvi_first, vvi_last);
    vvi_first = vvi_last;
  }

  // add bonds to virtual graph
  for (boost::tie(rei, rei_end) = bonds(rg); rei != rei_end; ++rei) {
    if (weight(*rei, rg)) {
      site_descriptor rs = source(*rei, rg);
      site_descriptor rt = target(*rei, rg);
      site_iterator vvsi, vvsi_end;
      for (boost::tie(vvsi, vvsi_end) = vm.virtual_sites(rg, rs);
           vvsi != vvsi_end; ++vvsi) {
        site_iterator vvti, vvti_end;
        for (boost::tie(vvti, vvti_end) = vm.virtual_sites(rg, rt);
             vvti != vvti_end; ++vvti) {
          bond_descriptor ved = add_edge(*vvsi, *vvti, vg).first;
          boost::put(bond_index_t(), vg, ved, num_bonds(vg) - 1);
          copy_property(bond_type_t(), rg, *rei, vg, ved);
          copy_property(boundary_crossing_t(), rg, *rei, vg, ved);
          copy_property(bond_vector_t(), rg, *rei, vg, ved);
          copy_property(bond_vector_relative_t(), rg, *rei, vg, ved);
        }
      }
    }
  }

  // add `in-real-site' bonds to virtual graph
  int dim = alps::get_or_default(dimension_t(), rg, int(0));
  alps::coordinate_type vec(dim, 0);
  for (boost::tie(rvi, rvi_end) = sites(rg); rvi != rvi_end; ++rvi) {
    if (weight(*rvi, rg)) {
      int t = vm.v2bond_type_offset() + boost::get(site_type_t(), rg, *rvi);
      site_iterator vvsi, vvsi_end;
      for (boost::tie(vvsi, vvsi_end) = vm.virtual_sites(rg, *rvi);
           vvsi != vvsi_end; ++vvsi) {
        for (site_iterator vvti = boost::next(vvsi); vvti != vvsi_end;
             ++vvti) {
          bond_descriptor ved = add_edge(*vvsi, *vvti, vg).first;
          boost::put(bond_index_t(), vg, ved, num_bonds(vg) - 1);
          boost::put(bond_type_t(), vg, ved, t);
          boost::put(boundary_crossing_t(), vg, ved, alps::boundary_crossing());
          boost::put(bond_vector_t(), vg, ved, vec);
          boost::put(bond_vector_relative_t(), vg, ved, vec);
        }
      }
    }
  }

  // setup bond and v2bond mapping
  bond_iterator vei_first = bonds(vg).first;
  bond_iterator vei_last = vei_first;
  for (boost::tie(rei, rei_end) = bonds(rg); rei != rei_end; ++rei) {
    if (weight(*rei, rg)) {
      vei_last += model.site(source(*rei, rg), rg).s().get_twice() *
        model.site(target(*rei, rg), rg).s().get_twice();
    }
    vm.add_bonds(rg, *rei, vei_first, vei_last);
    vei_first = vei_last;
  }

  for (boost::tie(rvi, rvi_end) = sites(rg); rvi != rvi_end; ++rvi) {
    if (weight(*rvi, rg)) {
      vei_last += model.site(*rvi, rg).s().get_twice() *
        (model.site(*rvi, rg).s().get_twice() - 1) / 2;
    }
    vm.add_v2bonds(rg, *rvi, vei_first, vei_last);
    vei_first = vei_last;
  }
}


//
// class template virtual_lattice
//

template<class GRAPH>
class virtual_lattice {
public:
  typedef GRAPH                       graph_type;
  typedef virtual_mapping<graph_type> mapping_type;

  typedef typename boost::graph_traits<graph_type>::vertex_descriptor
    vertex_descriptor;
  typedef typename boost::graph_traits<graph_type>::edge_descriptor
    edge_descriptor;
  typedef typename boost::graph_traits<graph_type>::adjacency_iterator
     adjacency_iterator;
  typedef typename boost::graph_traits<graph_type>::out_edge_iterator
     out_edge_iterator;
  typedef typename boost::graph_traits<graph_type>::in_edge_iterator
     in_edge_iterator;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::edge_iterator
    edge_iterator;

  typedef typename boost::graph_traits<graph_type>::directed_category
    directed_category;
  typedef typename boost::graph_traits<graph_type>::edge_parallel_category
    edge_parallel_category;
  typedef typename boost::graph_traits<graph_type>::traversal_category
    traversal_category;

  typedef typename boost::graph_traits<graph_type>::vertices_size_type
    vertices_size_type;
  typedef typename boost::graph_traits<graph_type>::edges_size_type
    edges_size_type;
  typedef typename boost::graph_traits<graph_type>::degree_size_type
    degree_size_type;

  typedef vertex_descriptor site_descriptor;
  typedef edge_descriptor   bond_descriptor;
  typedef vertex_iterator   site_iterator;
  typedef edge_iterator     bond_iterator;

  virtual_lattice() {}
  template<class RG, class M, class W>
  virtual_lattice(const RG& rg, const M& model, const W& weight)
  { initialize(rg, model, weight); }

  template<class RG>
  std::pair<site_iterator, site_iterator>
  virtual_sites(const RG& rg,
    const typename alps::graph_traits<RG>::site_descriptor& rv) const
  { return mapping_.virtual_sites(rg, rv); }
  template<class RG>
  std::pair<bond_iterator, bond_iterator>
  virtual_bonds(const RG& rg,
    const typename alps::graph_traits<RG>::bond_descriptor& re) const
  { return mapping_.virtual_bonds(rg, re); }

  void clear() { graph_.clear(); mapping_.clear(); }
  template<class RG, class M, class W>
  void initialize(const RG& rg, const M& model, const W& weight)
  { generate_virtual_lattice(graph_, mapping_, rg, model, weight); }

  const graph_type& graph() const { return graph_; }
  const mapping_type& mapping() const { return mapping_; }

  template<class RG>
  void print_mapping(std::ostream& os, const RG& rg) const {
    mapping_.output(os, rg, graph_);
  }

private:
  graph_type graph_;
  mapping_type mapping_;
};

}

namespace boost {

template<class RG>
std::pair<typename looper::virtual_lattice<RG>::vertex_iterator,
          typename looper::virtual_lattice<RG>::vertex_iterator>
vertices(const looper::virtual_lattice<RG>& vl)
{ return vertices(vl.graph()); }

template<class RG>
typename looper::virtual_lattice<RG>::vertices_size_type
num_vertices(const looper::virtual_lattice<RG>& vl)
{ return num_vertices(vl.graph()); }

template<class RG>
std::pair<typename looper::virtual_lattice<RG>::edge_iterator,
          typename looper::virtual_lattice<RG>::edge_iterator>
edges(const looper::virtual_lattice<RG>& vl)
{ return edges(vl.graph()); }

template<class RG>
typename looper::virtual_lattice<RG>::edges_size_type
num_edges(const looper::virtual_lattice<RG>& vl)
{ return num_edges(vl.graph()); }

} // end namespace boost

namespace looper {

template<class RG>
std::pair<typename virtual_lattice<RG>::site_iterator,
          typename virtual_lattice<RG>::site_iterator>
virtual_sites(const virtual_lattice<RG>& vl, const RG& rg,
  const typename alps::graph_traits<RG>::site_descriptor& rv)
{ return vl.virtual_sites(rg, rv); }

template<class RG>
std::pair<typename virtual_lattice<RG>::bond_iterator,
          typename virtual_lattice<RG>::bond_iterator>
virtual_bonds(const virtual_lattice<RG>& vl, const RG& rg,
  const typename alps::graph_traits<RG>::bond_descriptor& re)
{ return vl.virtual_bonds(rg, re); }

template<class RG>
void set_parity(virtual_lattice<RG>& vl)
{ alps::set_parity(vl.graph()); }

template<class RG>
int gauge(const virtual_lattice<RG>& vl,
  const typename virtual_lattice<RG>::site_descriptor& vd)
{ return gauge(vl.graph(), vd); }

template<class RG>
std::ostream& operator<<(std::ostream& os, const virtual_lattice<RG>& vl)
{
  os << vl.graph();
  return os;
}

} // end namespace looper

namespace alps {

template<class Graph>
struct parity_traits<looper::parity_t, Graph> {
  typedef int value_type;
  BOOST_STATIC_CONSTANT(value_type, white = 1);
  BOOST_STATIC_CONSTANT(value_type, black = -1);
  BOOST_STATIC_CONSTANT(value_type, undefined = 0);
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION

template<class Graph>
const typename parity_traits<looper::parity_t, Graph>::value_type
parity_traits<looper::parity_t, Graph>::white;

template<class Graph>
const typename parity_traits<looper::parity_t, Graph>::value_type
parity_traits<looper::parity_t, Graph>::black;

template<class Graph>
const typename parity_traits<looper::parity_t, Graph>::value_type
parity_traits<looper::parity_t, Graph>::undefined;

#endif

} // end namespace alps

#endif // LOOPER_LATTICE_H
