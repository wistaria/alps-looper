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

#include <looper/find_bridges.h>
#include <boost/graph/connected_components.hpp>
#include <boost/random.hpp>

struct bridge_visitor
{
  bridge_visitor(int& n) : nb(n) {}
  template<class Edge, class Graph>
  void operator()(Edge e, const Graph& g) const
  {
    std::cout << "  find a bridge (" << boost::source(e, g) << " -- "
              << boost::target(e,g) << ")\n";
    ++nb;
  }
  int& nb;
};

int main(int , char**)
{
  int nv;
  std::cin >> nv;

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::no_property,
    boost::property<boost::edge_index_t, std::size_t> >
    Graph;
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  Graph g(nv);

  // random number generator
  boost::variate_generator<boost::mt19937, boost::uniform_int<> >
    rng(boost::mt19937(8237), boost::uniform_int<>(0, nv-1));

  // make a connected graph
  std::vector<int> component(nv);
  do {
    int a = rng();
    int b = rng();
    boost::graph_traits<Graph>::edge_descriptor
      e = boost::add_edge(a,b,g).first;
    boost::put(boost::edge_index, g, e, boost::num_edges(g)-1);
    std::cout << "add edge: " << a << " -- " << b << std::endl;
  } while (boost::connected_components(g, &component[0]) > 1);

  typedef boost::default_color_type color_type;
  typedef boost::color_traits<color_type> color;
  std::vector<boost::default_color_type> vcolor_map(nv);
  std::vector<boost::default_color_type> ecolor_map(boost::num_edges(g));
  std::vector<boost::default_color_type> etype_map(boost::num_edges(g));
  std::vector<boost::graph_traits<Graph>::vertex_descriptor> in_postorder;

  boost::undirected_depth_first_visit(
    g,
    *boost::vertices(g).first,
    boost::make_dfs_visitor(
      std::make_pair(
        looper::stamp_edge(
          boost::make_iterator_property_map(
            etype_map.begin(),
            boost::get(boost::edge_index, g)),
          color::white(),
          boost::on_tree_edge()),
      std::make_pair(
        looper::stamp_edge(
          boost::make_iterator_property_map(
            etype_map.begin(),
            boost::get(boost::edge_index, g)),
          color::black(),
          boost::on_back_edge()),
        looper::write_vertex(
          std::back_inserter(in_postorder),
          boost::on_finish_vertex())))),
    boost::make_iterator_property_map(vcolor_map.begin(),
      boost::get(boost::vertex_index, g)),
    boost::make_iterator_property_map(ecolor_map.begin(),
      boost::get(boost::edge_index, g)));

  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
    if (etype_map[boost::get(boost::edge_index, g, *ei)] == color::white()) {
      std::cout << "tree edge (" << boost::source(*ei,g) << " -- "
                << boost::target(*ei,g) << ")\n";
    }
  }

  std::cout << "vertices in postorder: ";
  std::vector<boost::graph_traits<Graph>::vertex_descriptor>::const_iterator
    vi, vi_end;
  for (vi = in_postorder.begin(), vi_end = in_postorder.end();
       vi != vi_end; ++vi) {
    std::cout << ' ' << boost::get(boost::vertex_index, g, *vi);
  }
  std::cout << std::endl;

  int nb = 0;
  std::cout << "finding all the bridges in a naive method (O(N^2)):\n";
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
    if (boost::connected_components(
          looper::mask_edge(g, *ei),
          &component[0]) > 1) {
      std::cout << "  find a bridge (" << boost::source(*ei, g) << " -- "
                << boost::target(*ei,g) << ")\n";
      ++nb;
    }
  }
  std::cout << "  found " << nb << " bridge(s)\n";

  nb = 0;
  std::cout << "finding all the bridges by find_bridges() (O(N)):\n";
  looper::find_bridges(g, *boost::vertices(g).first, bridge_visitor(nb));
  std::cout << "  found " << nb << " bridge(s)\n";
}
