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

// $Id: virtual_graph.h 554 2003-11-12 02:36:24Z wistaria $

#ifndef LOOPER_VIRTUAL_GRAPH_H
#define LOOPER_VIRTUAL_GRAPH_H

#include <looper/graph.h>

#include <alps/model.h>
#include <utility>                 // std::pair, std::make_pair
#include <vector>                  // std::vector

namespace looper {

// class template vmapping
// for describing a mapping from a real site to virtual ones

template<class G>
class vmapping
{
public:
  typedef G                                           graph_type;
  typedef typename graph_type::vertex_iterator        vertex_iterator;
  typedef std::pair<vertex_iterator, vertex_iterator> range_type;

  vmapping() : map_(1, vertex_iterator()) {}
  
  int num_groups() const { return map_.size() - 1; }
  int num_virtual_vertices(int g) const {
    return *map_[g + 1] - *map_[g];
  }

  range_type virtual_vertices(int s) const {
    return std::make_pair(map_[s], map_[s + 1]);
  }

  void add(const vertex_iterator& first, const vertex_iterator& last) {
    map_[map_.size() - 1] = first;
    map_.push_back(last);
  }
  void clear() {
    map_.clear();
    map_.push_back(vertex_iterator());
  }

  friend bool operator==(const vmapping& lhs,
			 const vmapping& rhs) {
    return lhs.map_ == rhs.map_;
  }
  friend bool operator!=(const vmapping& lhs,
			 const vmapping& rhs) {
    return !operator==(lhs, rhs);
  }

  void output(std::ostream& os) const {
    os << "[[vmapping]]\n";
    os << "  number of groups = " << num_groups() << std::endl;
    os << "  mapping:\n";
    for (int g = 0; g < num_groups(); ++g) {
      os << "    " << g << " -> ";
      if (num_virtual_vertices(g) == 0) {
	os << "null\n";
      } else if (num_virtual_vertices(g) == 1) {
	os << *(map_[g]) << '\n';
      } else {
	os << "[" << *(map_[g]) << "," << *(map_[g+1]) - 1 << "]\n";
      }
    }
  }

  friend std::ostream& operator<<(std::ostream& os,
				  const vmapping& map) {
    map.output(os);
    return os;
  }
  
private:
  std::vector<vertex_iterator> map_;
};

  
// class template virtual_graph

template<class G>
struct virtual_graph
{
  typedef G           graph_type;
  typedef vmapping<G> mapping_type;

  graph_type   graph;
  mapping_type mapping;
  int num_real_vertices;
  int num_real_edges;
};


// function generate_virtual_graph

template<class RG, class MDL, class G>
inline void generate_virtual_graph(const RG& rg, const MDL& model,
				   virtual_graph<G>& vg)
{
  typedef RG                                      rgraph_type;
  typedef typename virtual_graph<G>::graph_type   vgraph_type;
  typedef typename virtual_graph<G>::mapping_type mapping_type;

  vg.graph.clear();
  vg.mapping.clear();
  vg.num_real_vertices = boost::num_vertices(rg);
  vg.num_real_edges = boost::num_edges(rg);

  // setup graph properties
  boost::get_property(vg.graph, graph_name_t())
    = "virtual graph of " + boost::get_property(rg, graph_name_t());
  boost::get_property(vg.graph, dimension_t())
    = boost::get_property(rg, dimension_t());

  // setup vertices
  typename rgraph_type::vertex_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = boost::vertices(rg); rvi != rvi_end; ++rvi) {
    int t = boost::get(vertex_type_t(), rg, *rvi);
    for (int i = 0; i < model.spin(t).get_twice(); ++i) {
      // add vertices to virtual graph
      typename vgraph_type::vertex_descriptor
	vvd = boost::add_vertex(vg.graph);

      // copy vertex properties
      alps::copy_property(vertex_type_t(), rg, *rvi, vg.graph, vvd);
      alps::copy_property(coordinate_t(), rg, *rvi, vg.graph, vvd);
      alps::copy_property(parity_t(), rg, *rvi, vg.graph, vvd);
    }
  }

  // setup mapping
  typename vgraph_type::vertex_iterator
    vvi_first = boost::vertices(vg.graph).first;
  typename vgraph_type::vertex_iterator vvi_last = vvi_first;
  for (boost::tie(rvi, rvi_end) = boost::vertices(rg); rvi != rvi_end; ++rvi) {
    int t = boost::get(vertex_type_t(), rg, *rvi);
    for (int i = 0; i < model.spin(t).get_twice(); ++i)
      ++vvi_last;
    vg.mapping.add(vvi_first, vvi_last);
    vvi_first = vvi_last;
  }

  // setup edges
  typename rgraph_type::edge_iterator rei, rei_end;
  for (boost::tie(rei, rei_end) = boost::edges(rg); rei != rei_end; ++rei) {
    typename rgraph_type::vertex_descriptor rs = boost::source(*rei, rg);
    typename rgraph_type::vertex_descriptor rt = boost::target(*rei, rg);
    typename vgraph_type::vertex_iterator vvsi, vvsi_end;
    for (boost::tie(vvsi, vvsi_end) = vg.mapping.virtual_vertices(rs);
	 vvsi != vvsi_end; ++vvsi) {
      typename vgraph_type::vertex_iterator vvti, vvti_end;
      for (boost::tie(vvti, vvti_end) = vg.mapping.virtual_vertices(rt);
	   vvti != vvti_end; ++vvti) {
	// add edges to virtual graph
	typename vgraph_type::edge_descriptor ved =
	  boost::add_edge(*vvsi, *vvti, vg.graph).first;
	  
	// setup edge properties
	boost::put(edge_index_t(), vg.graph, ved, boost::num_edges(vg.graph));
	alps::copy_property(edge_type_t(), rg, *rei, vg.graph, ved);
      }
    }
  }
}

namespace vg_detail {

template<class T>
struct vector_spin_wrapper
{
  vector_spin_wrapper(const std::vector<T>& v) : vec_(v) {}
  const T& spin(int i) const { return vec_[i]; }
  const std::vector<T>& vec_;
};

template<class T>
struct const_spin_wrapper
{
  const_spin_wrapper(const T& t) : t_(t) {}
  const T& spin(int) const { return t_; }
  const T& t_;
};

}

template<class G, class IntType>
inline void generate_virtual_graph(const G& rg,
				   alps::half_integer<IntType> s,
				   virtual_graph<G>& vg)
{
  generate_virtual_graph(rg,
    vg_detail::const_spin_wrapper<alps::half_integer<IntType> >(s), vg);
}

template<class G, class IntType>
inline void generate_virtual_graph(const G& rg,
				   std::vector<alps::half_integer<IntType> > v,
				   virtual_graph<G>& vg)
{
  generate_virtual_graph(rg,
    vg_detail::vector_spin_wrapper<alps::half_integer<IntType> >(v), vg);
}
  
} // end namespace looper

#endif // LOOPER_VIRTUAL_GRAPH_H
