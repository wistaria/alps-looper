/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: graph.h 561 2003-11-12 15:53:43Z wistaria $

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
                                              alps::detail::coordinate_type,
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
                                              alps::detail::coordinate_type,
                                boost::property<parity_t, int8_t,
                                boost::property<vertex_type_t, int > > >,
                              boost::property<edge_type_t, int,
                                boost::property<edge_index_t, int > >,
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

template<class D = unsigned int, class S = D, class E = std::vector<S> >
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

namespace detail {

template<class G, class D>
struct hypercubic_graph_helper
{
  typedef G                                                  graph_type;
  typedef D                                                  descriptor_type;
  typedef typename descriptor_type::dimension_type           dimension_type;
  typedef typename descriptor_type::size_type                index_type;
  typedef typename descriptor_type::extent_type              base_type;
  typedef base_type                                          position_type;
  typedef typename graph_traits<graph_type>::coordinate_type coordinate_type;

  template<class E>
  static position_type index2pos(const E& ext, const base_type& base,
                                 index_type i) {
    position_type pos(ext.size());
    for (dimension_type d = 0; d < ext.size(); ++d) {
      pos[d] = (i / base[d]) % ext[d];
    }
    return pos;
  }

  static index_type pos2index(dimension_type dim, const base_type& base,
                              const position_type& pos) {
    index_type i = 0;
    for (dimension_type d = 0; d < dim; ++d) i += base[d] * pos[d];
    return i;
  }

  template<class E>
  static std::size_t neighbor(const E& ext, const position_type& base,
                              dimension_type s, dimension_type d,
                              typename position_type::value_type p) {
    position_type pos = index2pos(ext, base, s);
    pos[d] = (pos[d] + ext[d] + p) % ext[d];
    return pos2index(ext.size(), base, pos);
  }

  static coordinate_type pos2coord(const position_type& pos) {
    coordinate_type coord(pos.size());
    for (dimension_type i = 0; i < pos.size(); ++i) coord[i] = pos[i];
    return coord;
  }
};

} // end namespace detail

template<class D, class S, class E, class G>
void generate_graph(const hypercubic_graph_generator<D, S, E>& desc,
                    G& g)
{
  typedef G graph_type;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::vertex_descriptor
    vertex_descriptor;
  typedef typename boost::graph_traits<graph_type>::edge_descriptor
    edge_descriptor;
  typedef typename graph_traits<graph_type>::coordinate_type coordinate_type;

  typedef hypercubic_graph_generator<D, S, E>                descriptor_type;
  typedef detail::hypercubic_graph_helper<graph_type, descriptor_type>
                                                             helper_type;

  typedef typename descriptor_type::dimension_type           dimension_type;
  typedef typename descriptor_type::size_type                size_type;
  typedef typename helper_type::position_type                position_type;

  if (desc.dimension() == 0) 
    boost::throw_exception(std::runtime_error("invalid dimension"));
  for (dimension_type d = 0; d < desc.dimension(); ++d) {
    if (desc.length(d) == 0)
      boost::throw_exception(std::runtime_error("invalid extent"));
  }

  g.clear();
  boost::get_property(g, graph_name_t()) = "simple hypercubic graph";
  boost::get_property(g, dimension_t()) = desc.dimension();

  size_type size = 1;
  position_type base(desc.dimension());
  for (dimension_type d = 0; d < desc.dimension(); ++d) {
    base[d] = size;
    size *= desc.length(d);
  }
  
  // setup vertices
  for (size_type i = 0; i < size; ++i) {
    vertex_descriptor vd = boost::add_vertex(g);
    coordinate_type coord =
      helper_type::pos2coord(helper_type::index2pos(desc.extent(), base, i));
    boost::put(vertex_type_t(), g, vd, 0);
    boost::put(coordinate_t(), g, vd, coord);
  }
  
  // setup edges
  for (size_type i = 0; i < size; ++i) {
    vertex_iterator vsi = boost::vertices(g).first + i;
    for (dimension_type d = 0; d < desc.dimension(); ++d) {
      if (desc.length(d) == 1) {
        // nothing to do
      } else {
        // avoid double bond
        if (desc.length(d) > 2 ||
            helper_type::index2pos(desc.extent(), base, i)[d] + 1
              != desc.length(d)) {
          vertex_iterator vti = boost::vertices(g).first +
            helper_type::neighbor(desc.extent(), base, i, d, 1);
          edge_descriptor ed = boost::add_edge(*vsi, *vti, g).first;
          boost::put(edge_index_t(), g, ed, boost::num_edges(g));
          boost::put(edge_type_t(), g, ed, d);
        }
      }
    }
  }
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
