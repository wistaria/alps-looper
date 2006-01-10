/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004 by Synge Todo <wistaria@comp-phys.org>,
*                       Marc Locher <Marc@itp.phys.ethz.ch>
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

// Algorithm for finding all the bridges on a connected graph

// Reference:
//   R. E. Tarjan, A Note on Finding the Bridges on a Graph,
//   Information Processing Letters 2 160-161 (1974).

#ifndef LOOPER_FINDBRIDGE_H
#define LOOPER_FINDBRIDGE_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/undirected_dfs.hpp>

namespace looper {

// vertex_writer & write_vertex

template<class OutputIterator, class Tag>
struct vertex_writer
  : public boost::base_visitor<vertex_writer<OutputIterator, Tag> >
{
  typedef Tag event_filter;
  vertex_writer(OutputIterator out) : m_out(out) {}
  template <class Vertex, class Graph>
  void operator()(Vertex v, const Graph&) { *m_out++ = v; }
  OutputIterator m_out;
};
template<class OutputIterator, class Tag>
inline vertex_writer<OutputIterator, Tag>
write_vertex(OutputIterator out, Tag)
{
  return vertex_writer<OutputIterator, Tag>(out);
}


// edge_stamper & stamp_edge

template<class PMap, class PValue, class Tag>
struct edge_stamper
  : public boost::base_visitor<edge_stamper<PMap, PValue, Tag> >
{
  typedef Tag event_filter;
  edge_stamper(PMap pm, const PValue& pv) : pmap(pm), pval(pv) {}
  template <class Edge, class Graph>
  void operator()(Edge e, const Graph&) { boost::put(pmap, e, pval); }
  PMap pmap;
  PValue pval;
};
template<class PMap, class PValue, class Tag>
inline edge_stamper<PMap, PValue, Tag>
stamp_edge(const PMap& pm, const PValue& pv, Tag)
{
  return edge_stamper<PMap, PValue, Tag>(pm, pv);
}


// edge_mask & mask_edge

template<typename Graph>
struct edge_mask
{
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
  edge_mask() : ex() {}
  edge_mask(Edge e) : ex(e) {}
  bool operator()(Edge e) const { return (e != ex); }
  Edge ex;
};
template<typename Graph>
inline boost::filtered_graph<Graph, edge_mask<Graph> >
mask_edge(
  const Graph& g,
  typename boost::graph_traits<Graph>::edge_descriptor e)
{
  return boost::filtered_graph<Graph, edge_mask<Graph> >(
    g,
    edge_mask<Graph>(e));
}


// algorithm find_bridges

template<typename Graph, typename BridgeVisitor>
inline void find_bridges(
  const Graph& g,
  typename boost::graph_traits<Graph>::vertex_descriptor start_vertex,
  BridgeVisitor vis)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  int nv = boost::num_vertices(g);
  int ne = boost::num_edges(g);

  typedef boost::default_color_type color_type;
  typedef boost::color_traits<color_type> color;
  typedef std::vector<color_type> color_map_type;
  color_map_type vcolor_map(nv);
  color_map_type ecolor_map(ne);
  color_map_type etype_map(ne);
  std::vector<Vertex> in_postorder;

  boost::undirected_depth_first_visit(
    g,
    start_vertex,
    boost::make_dfs_visitor(
      std::make_pair(
        boost::write_property(
          boost::get(boost::vertex_index, g),
          std::back_inserter(in_postorder),
          boost::on_finish_vertex()),
      std::make_pair(
        stamp_edge(
          boost::make_iterator_property_map(
            etype_map.begin(),
            boost::get(boost::edge_index, g)),
          color::white(),
          boost::on_tree_edge()),
        stamp_edge(
          boost::make_iterator_property_map(
            etype_map.begin(),
            boost::get(boost::edge_index, g)),
          color::black(),
          boost::on_back_edge())))),
    boost::make_iterator_property_map(vcolor_map.begin(),
      boost::get(boost::vertex_index, g)),
    boost::make_iterator_property_map(ecolor_map.begin(),
      boost::get(boost::edge_index, g)));

  std::vector<int> ND(nv);
  std::vector<int> L(nv);
  std::vector<int> H(nv);
  std::vector<int> reverse_postorder(nv);

  typename std::vector<Vertex>::const_iterator vi, vi_end;
  int v;
  for (vi = in_postorder.begin(), vi_end = in_postorder.end(), v = 0;
       vi != vi_end; ++vi, ++v) reverse_postorder[*vi] = v;

  for (vi = in_postorder.begin(), vi_end = in_postorder.end(), v = 0;
       vi != vi_end; ++vi, ++v) {
    ND[v] = 1;
    L[v] = nv;
    H[v] = 0;
    typename boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::out_edges(*vi, g);
         ei != ei_end; ++ei) {
      int w = reverse_postorder[boost::opposite(*ei, *vi, g)];
      if (etype_map[boost::get(boost::edge_index,g,*ei)] == color::white()) {
        // tree edge
        if (w < v) {
          // w is a son of v
          ND[v] += ND[w];
          L[v] = std::min(L[v], L[w]);
          H[v] = std::max(H[v], H[w]);
          // check if this tree edge is bridge or not
          if ((H[w] <= w) && (L[w] > w - ND[w])) vis(*ei, g);
        } else {
          // w is the father of v, skip
        }
      } else {
        if (w != v) {
          // back edge
          L[v] = std::min(L[v], w);
          H[v] = std::max(H[v], w);
        }
      }
    }
    L[v] = std::min(L[v], v - ND[v] + 1);
    H[v] = std::max(H[v], v);
    // std::cout << v << ' ' << in_postorder[v] << ' ' << ND[v] << ' '
    //           << L[v] << ' ' << H[v] << std::endl;
  }
}


// find_bridges with a skip edge

template<typename Graph, typename BridgeVisitor>
inline void find_bridges(
  const Graph& g,
  typename boost::graph_traits<Graph>::vertex_descriptor start_vertex,
  BridgeVisitor& vis,
  typename boost::graph_traits<Graph>::edge_descriptor skip_edge)
{
  find_bridges(mask_edge(g, skip_edge), start_vertex, vis);
}

} // end namespace looper

#endif // LOOPER_FINDBRIDGE_H
