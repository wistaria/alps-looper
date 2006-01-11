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

#include <looper/find_bridge.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <cmath>
#include <boost/random.hpp>

template<class Graph>
struct my_visitor : public boost::dfs_visitor<> {
  typedef boost::dfs_visitor<> BASE_;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

  my_visitor() : BASE_() {}

  void initialize_vertex(Vertex v, const Graph&) {
    std::cout << "initialize_vertex: " << v << std::endl;
  }
  void discover_vertex(Vertex v, const Graph&) {
    std::cout << "discover_vertex: " << v << std::endl;
  }
  void finish_vertex(Vertex v, const Graph&) {
    std::cout << "finish_vertex: " << v << std::endl;
  }
  void examine_edge(Edge e, const Graph& g) {
    std::cout << "examine_edge: " << source(e, g) << " -- "
              << target(e, g) << std::endl;
  }
  void back_edge(Edge e, const Graph& g) {
    std::cout << "back_edge: " << source(e, g) << " -- "
              << target(e, g) << std::endl;
  }
  void tree_edge(Edge e, const Graph& g) {
    std::cout << "tree_edge: " << source(e, g) << " -- "
              << target(e, g) << std::endl;
  }
};

int main(int , char**)
{
  int nv;
  std::cin >> nv;

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property,
    boost::property<boost::edge_index_t, std::size_t> >
    Graph;
  typedef boost::graph_traits<Graph>::edge_descriptor Edge;
  Graph g(nv);

  // random number generator
  boost::variate_generator<boost::mt19937, boost::uniform_int<> >
    rng(boost::mt19937(8237), boost::uniform_int<>(0, nv-1));

  int a0 = 0;
  int b0 = 1;
  Edge e = add_edge(a0,b0,g).first;
  put(boost::edge_index, g, e, num_edges(g)-1);
  std::cout << "add edge: " << a0 << " -- " << b0 << std::endl;
  std::vector<int> component(nv);
  do {
    int a = rng();
    int b = rng();
    if ((a != a0 || b != b0) && (a != b0 || b != a0)) {
      Edge e = add_edge(a,b,g).first;
      put(boost::edge_index, g, e, num_edges(g)-1);
      std::cout << "add edge: " << a << " -- " << b << std::endl;
    }
  } while (connected_components(g, &component[0]) > 1);

  std::vector<boost::default_color_type> vcolor_map(num_vertices(g));
  std::vector<boost::default_color_type> ecolor_map(num_edges(g));

  std::cout << "DFS on original graph\n";
  undirected_dfs(
    g,
    my_visitor<Graph>(),
    boost::make_iterator_property_map(vcolor_map.begin(),
      get(boost::vertex_index, g)),
    boost::make_iterator_property_map(ecolor_map.begin(),
      get(boost::edge_index, g)));

  std::cout << "DFS on filtered graph (edge 0 -- 1 is removed)\n";
  undirected_dfs(
    looper::mask_edge(g, *edges(g).first),
    my_visitor<boost::filtered_graph<Graph, looper::edge_mask<Graph> > >(),
    boost::make_iterator_property_map(vcolor_map.begin(),
      get(boost::vertex_index,
        looper::mask_edge(g, *edges(g).first))),
    boost::make_iterator_property_map(ecolor_map.begin(),
      get(boost::edge_index,
        looper::mask_edge(g, *edges(g).first))));

  return 0;
}
