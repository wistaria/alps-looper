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

#ifndef LOOPER_GRAPH_H__
#define LOOPER_GRAPH_H__

#include <alps/lattice.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>
#include <vector>

namespace looper {

typedef alps::graph_name_t    graph_name_t;
typedef alps::dimension_t     dimension_t;
typedef boost::vertex_index_t vertex_index_t;
typedef alps::vertex_type_t   vertex_type_t;
typedef alps::coordinate_t    coordinate_t;
typedef alps::parity_t        parity_t;
typedef boost::edge_index_t   edge_index_t;
typedef alps::edge_type_t     edge_type_t;

// NOTE: We use adjacency_list with EdgeListS=vecS, which will be less
// efficient than that with EdgeListS=listS for edge insertion.
// However, the former allows us random access to edges, which is
// required for QMC in SSE representation.

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS,
                              boost::property<coordinate_t,
                                              alps::coordinate_type,
                                boost::property<vertex_type_t, int > >,
                              boost::property<edge_type_t, int,
                                boost::property<edge_index_t, int > >,
                              boost::property<dimension_t, std::size_t,
                                boost::property<graph_name_t, std::string > >,
                              boost::vecS> graph_type;

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS,
                              boost::property<coordinate_t,
                                              alps::coordinate_type,
                                boost::property<parity_t, int,
                                boost::property<vertex_type_t, int> > >,
                              boost::property<edge_type_t, int,
                                boost::property<edge_index_t, int> >,
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

  // unit cell
  alps::lattice_traits<lattice_type>::unit_cell_type uc(dim);
  uc.add_vertex(0, alps::lattice_traits<lattice_type>::unit_cell_type::
                coordinate_type(dim, 0.0));
  alps::lattice_traits<lattice_type>::offset_type so(dim, 0), to(dim, 0);
  for (std::size_t d = 0; d < dim; ++d) {
    to[d] = 1; uc.add_edge(d, 1, so, 1, to); to[d] = 0;
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
inline int
site_type(typename boost::graph_traits<
            boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
            vertex_descriptor vd,
          const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  return boost::get(alps::site_type_t(), g, vd);
}

template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
inline int
bond_type(typename boost::graph_traits<
            boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6> >::
            edge_descriptor ed,
          const boost::adjacency_list<T0, T1, T2, T3, T4, T5, T6>& g)
{
  return boost::get(alps::bond_type_t(), g, ed);
}

} // end namespace looper

#endif // LOOPER_GRAPH_H
