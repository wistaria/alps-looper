/**************************************************************************** 
*
* alps/looper: multi-cluster quantum Monte Carlo algorithm for spin systems
*              in path-integral and SSE representations
*
* $Id: virtualgraph.h 398 2003-10-09 10:33:05Z wistaria $
*
* Copyright (C) 2001-2003 by Synge Todo <wistaria@comp-phys.org>,
*
* Permission is hereby granted, free of charge, to any person or organization 
* obtaining a copy of the software covered by this license (the "Software") 
* to use, reproduce, display, distribute, execute, and transmit the Software, 
* and to prepare derivative works of the Software, and to permit others
* to do so for non-commerical academic use, all subject to the following:
*
* The copyright notice in the Software and this entire statement, including 
* the above license grant, this restriction and the following disclaimer, 
* must be included in all copies of the Software, in whole or in part, and 
* all derivative works of the Software, unless such copies or derivative 
* works are solely in the form of machine-executable object code generated by 
* a source language processor.
*
* In any scientific publication based in part or wholly on the Software, the
* use of the Software has to be acknowledged and the publications quoted
* on the web page http://www.alps.org/license/ have to be referenced.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
**************************************************************************/

#ifndef LOOPER_VIRTUAL_GRAPH_H
#define LOOPER_VIRTUAL_GRAPH_H

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

  
// class template virtual_graph_type
// for choosing proper type for virtual graph and mapping

template<class G>
struct virtual_graph_type
{
  typedef G                    graph_type;
  typedef vmapping<graph_type> mapping_type;
};


// function generate_virtual_graph
// for generating virtual graph and mapping from real graph

template<class G, class IntType>
inline void generate_virtual_graph(const G& rg,
  const std::vector<alps::half_integer<IntType> >& spins,
  typename virtual_graph_type<G>::graph_type& vg,
  typename virtual_graph_type<G>::mapping_type& mp)
{
  typedef G                                            rgraph_type;
  typedef alps::half_integer<IntType>                  spin_type;
  typedef typename virtual_graph_type<G>::graph_type   vgraph_type;
  typedef typename virtual_graph_type<G>::mapping_type mapping_type;

  vg.clear();
  mp.clear();

  // setup graph properties
  boost::get_property(vg, graph_name_t())
    = "virtual graph of " + boost::get_property(rg, graph_name_t());
  boost::get_property(vg, dimension_t())
    = boost::get_property(rg, dimension_t());

  // setup vertices
  typename rgraph_type::vertex_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = boost::vertices(rg); rvi != rvi_end; ++rvi) {
    int t = boost::get(vertex_type_t(), rg, *rvi);
    for (int i = 0; i < spins.at(t).get_twice(); ++i) {
      // add vertices to virtual graph
      typename vgraph_type::vertex_descriptor vvd = boost::add_vertex(vg);

      // copy vertex properties
      alps::copy_property(vertex_type_t(), rg, *rvi, vg, vvd);
      alps::copy_property(coordinate_t(), rg, *rvi, vg, vvd);
      alps::copy_property(parity_t(), rg, *rvi, vg, vvd);
    }
  }

  // setup mapping
  typename vgraph_type::vertex_iterator vvi_first = boost::vertices(vg).first;
  typename vgraph_type::vertex_iterator vvi_last = vvi_first;
  for (boost::tie(rvi, rvi_end) = boost::vertices(rg); rvi != rvi_end; ++rvi) {
    int t = boost::get(vertex_type_t(), rg, *rvi);
    for (int i = 0; i < spins.at(t).get_twice(); ++i)
      ++vvi_last;
    mp.add(vvi_first, vvi_last);
    vvi_first = vvi_last;
  }

  // setup edges
  typename rgraph_type::edge_iterator rei, rei_end;
  for (boost::tie(rei, rei_end) = boost::edges(rg); rei != rei_end; ++rei) {
    typename rgraph_type::vertex_descriptor rs = boost::source(*rei, rg);
    typename rgraph_type::vertex_descriptor rt = boost::target(*rei, rg);
    typename vgraph_type::vertex_iterator vvsi, vvsi_end;
    for (boost::tie(vvsi, vvsi_end) = mp.virtual_vertices(rs);
	 vvsi != vvsi_end; ++vvsi) {
      typename vgraph_type::vertex_iterator vvti, vvti_end;
      for (boost::tie(vvti, vvti_end) = mp.virtual_vertices(rt);
	   vvti != vvti_end; ++vvti) {
	// add edges to virtual graph
	typename vgraph_type::edge_descriptor ved =
	  boost::add_edge(*vvsi, *vvti, vg).first;
	  
	// setup edge properties
	boost::put(edge_index_t(), vg, ved, boost::num_edges(vg));
	alps::copy_property(edge_type_t(), rg, *rei, vg, ved);
      }
    }
  }
}

template<class G, class IntType>
inline void generate_virtual_graph(const G& rg, alps::half_integer<IntType> s,
  typename virtual_graph_type<G>::graph_type& vg,
  typename virtual_graph_type<G>::mapping_type& mp)
{
  int maxtype = 0;
  typename G::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(rg); vi != vi_end; ++vi)
    maxtype = std::max(maxtype, boost::get(vertex_type_t(), rg, *vi));
  std::vector<alps::half_integer<IntType> > spins(maxtype + 1, s);
  generate_virtual_graph(rg, spins, vg, mp);
}
  
} // end namespace looper

#endif // LOOPER_VIRTUAL_GRAPH_H
