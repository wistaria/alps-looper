/**************************************************************************** 
*
* alps/looper: multi-cluster quantum Monte Carlo algorithm for spin systems
*              in path-integral and SSE representations
*
* $Id: sse.h 486 2003-10-30 05:01:15Z wistaria $
*
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>,
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

#ifndef LOOPER_SSE_H
#define LOOPER_SSE_H

#include <looper/permutation.h>
#include <looper/union_find.h>
#include <looper/virtual_graph.h>
#include <looper/weight.h>
#include <boost/throw_exception.hpp>
#include <cmath>
#include <stdexcept>

namespace looper {

typedef qmc_node sse_node;


template<class G, class M, class W = default_weight> struct sse;

template<class G, class M, class W>
struct sse<virtual_graph<G>, M, W>
{
  typedef virtual_graph<G>                      vg_type;
  typedef typename virtual_graph<G>::graph_type graph_type;
  typedef M                                     model_type;
  typedef W                                     weight_type;

  typedef typename boost::graph_traits<graph_type>::edge_iterator
                                                    edge_iterator;
  typedef typename boost::graph_traits<graph_type>::edge_descriptor
                                                    edge_descriptor;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
                                                    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::vertex_descriptor
                                                    vertex_descriptor;

  struct parameter_type
  {
    typedef virtual_graph<G> vg_type;
    typedef typename vg_type::graph_type graph_type;
    typedef typename vg_type::mapping_type mapping_type;
    typedef M model_type;
    typedef W weight_type;

    parameter_type(const vg_type& vg, const model_type& m, double b)
      : virtual_graph(vg), graph(vg.graph), mapping(vg.mapping), model(m),
	beta(b), chooser(vg, model) {}

    const vg_type&                  virtual_graph;
    const graph_type&               graph;
    const mapping_type&             mapping;
    const model_type&               model;
    const double                    beta;
    const bond_chooser<weight_type> chooser;
  };

  struct config_type
  {
    typedef std::vector<sse_node>            os_type;
    typedef typename os_type::iterator       iterator;
    typedef typename os_type::const_iterator const_iterator;
    typedef typename os_type::value_type     node_type;

    os_type bottom;
    os_type top;
    os_type os;
    unsigned int num_operators;

    unsigned int num_loops0;
    unsigned int num_loops;
  };

  typedef typename config_type::iterator       operator_iterator;
  typedef typename config_type::const_iterator const_operator_iterator;
  
  //
  // update functions
  //
  
  // initialize
  static void initialize(config_type& config, const vg_type& vg, int ni = 16)
  {
    config.bottom.clear();
    config.bottom.resize(boost::num_vertices(vg.graph));
    config.top.clear();
    config.top.resize(boost::num_vertices(vg.graph));

    vertex_iterator vi_end = boost::vertices(vg.graph).second;
    for (vertex_iterator vi = boost::vertices(vg.graph).first;
	 vi != vi_end; ++vi) {
      // all up
      config.bottom[*vi].conf() = 0;
      config.top[*vi].conf() = 0;
    }
    config.os.clear();
    config.os.resize(ni);
    operator_iterator oi_end = config.os.end();
    for (operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
      oi->clear_graph();
      oi->set_to_identity();
    }
    config.num_operators = 0;
  }
  static void initialize(config_type& config, const parameter_type& p,
			 int ni = 16)
  { initialize(config, p.virtual_graph, ni); }

  template<class RNG>
  static void generate_loops(config_type& config, const vg_type& vg,
			     const model_type& /* model */, double beta,
			     const bond_chooser<weight_type>& bc,
			     RNG& uniform_01)
  {
    //
    // diagonal update & labeling
    //

    // copy spin configurations at the bottom
    std::vector<int> curr_conf(boost::num_vertices(vg.graph));
    {
      vertex_iterator vi, vi_end;
      std::vector<int>::iterator itr = curr_conf.begin();
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph); 
	   vi != vi_end; ++vi, ++itr) *itr = config.bottom[*vi].conf();
    }
    
    // scan over operators
    {
      operator_iterator oi_end = config.os.end();
      for (operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
	if (oi->is_identity()) {
	  // identity operator
	  int b = bc.choose(uniform_01);
	  edge_iterator ei = boost::edges(vg.graph).first + b;
	  if (uniform_01() <
	      bc.global_weight() * beta *
	      bc.weight(b).p_accept(curr_conf[boost::source(*ei, vg.graph)],
				    curr_conf[boost::target(*ei, vg.graph)]) /
	      double(config.os.size() - config.num_operators)) {
	    // insert diagonal operator
	    oi->identity_to_diagonal();
	    oi->set_bond(b);
	    ++config.num_operators;
	    
	    oi->set_new((curr_conf[boost::source(*ei, vg.graph)] ^
			 curr_conf[boost::target(*ei, vg.graph)]),
			(uniform_01() < bc.weight(b).p_freeze()));
	  } else { /* nothing to be done */ }
	} else if (oi->is_diagonal()) {
	  // diagonal operator
	  int b = oi->bond();
	  edge_iterator ei = boost::edges(vg.graph).first + b;
	  if (uniform_01() <
	      double(config.os.size() - config.num_operators + 1) / 
	      (bc.global_weight() * beta *
	       bc.weight(b).
	         p_accept(curr_conf[boost::source(*ei, vg.graph)],
			  curr_conf[boost::target(*ei, vg.graph)]))) {
	    // remove diagonal operator
	    oi->diagonal_to_identity();
	    --config.num_operators;
	  } else {
	    // remain as diagonal operator
	    oi->set_old((curr_conf[boost::source(*ei, vg.graph)] ^
			 curr_conf[boost::target(*ei, vg.graph)]));
	  }
	} else {
	  // off-diagonal operator
	  int b = oi->bond();
	  edge_iterator ei = boost::edges(vg.graph).first + oi->bond();
	  curr_conf[boost::source(*ei, vg.graph)] ^= 1;
	  curr_conf[boost::target(*ei, vg.graph)] ^= 1;
	  oi->set_old(uniform_01() < bc.weight(b).p_reflect());
	}
      }
    }
#ifndef NDEBUG
    // check consistency of spin configuration
    {
      vertex_iterator vi, vi_end;
      std::vector<int>::iterator itr = curr_conf.begin();
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph); 
	   vi != vi_end; ++vi, ++itr) assert(*itr == config.top[*vi].conf());
    }
#endif

    //
    // cluster identification by using union-find algorithm
    //

    {
      std::vector<sse_node::segment_type *>
	curr_ptr(boost::num_vertices(vg.graph));
      {
	vertex_iterator vi, vi_end;
	std::vector<sse_node::segment_type *>::iterator pi = curr_ptr.begin();
	for (boost::tie(vi, vi_end) = boost::vertices(vg.graph); 
	     vi != vi_end; ++vi, ++pi)
	  *pi = &(config.bottom[*vi].loop_segment(0));
      }

      // scan over operators
      {
	operator_iterator oi_end = config.os.end();
	for (operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
	  if (!oi->is_identity()) {
	    edge_iterator ei = boost::edges(vg.graph).first + oi->bond();
	    union_find::unify(*curr_ptr[boost::source(*ei, vg.graph)],
			      segment_d(oi, 0));
	    union_find::unify(*curr_ptr[boost::target(*ei, vg.graph)],
			      segment_d(oi, 1));
	    if (oi->is_frozen())
	      union_find::unify(oi->loop_segment(0), oi->loop_segment(1));
	    curr_ptr[boost::source(*ei, vg.graph)] = &segment_u(oi, 0);
	    curr_ptr[boost::target(*ei, vg.graph)] = &segment_u(oi, 1);
	  }
	}
      }

      // connect to top
      {
	vertex_iterator vi, vi_end;
	for (boost::tie(vi, vi_end) = boost::vertices(vg.graph); 
	     vi != vi_end; ++vi)
	  union_find::unify(*curr_ptr[*vi], config.top[*vi].loop_segment(0));
      }
    }

    // connect bottom and top with random permutation
    {
      std::vector<int> r, c0, c1;
      for (int i = 0; i < vg.mapping.num_groups(); ++i) {
	int s2 = vg.mapping.num_virtual_vertices(i);
	int offset = *(vg.mapping.virtual_vertices(i).first);
	if (s2 == 1) {
	  // S=1/2: just connect top and bottom
	  union_find::unify(config.bottom[offset].loop_segment(0),
			    config.top   [offset].loop_segment(0));
	} else {
	  // S>=1
	  r.resize(s2);
	  c0.resize(s2);
	  c1.resize(s2);
	  vertex_iterator vi_end = vg.mapping.virtual_vertices(i).second;
	  for (vertex_iterator vi = vg.mapping.virtual_vertices(i).first;
	       vi != vi_end; ++vi) {
	    r[*vi - offset] = *vi - offset;
	    c0[*vi - offset] = config.bottom[*vi].conf();
	    c1[*vi - offset] = config.top   [*vi].conf();
	  }
	  restricted_random_shuffle(r.begin(), r.end(),
				    c0.begin(), c0.end(),
				    c1.begin(), c1.end(),
				    uniform_01);
	  for (int j = 0; j < s2; ++j) {
	    union_find::unify(
			      config.bottom[offset +   j ].loop_segment(0),
			      config.top   [offset + r[j]].loop_segment(0));
	  }
	}
      }
    }

    // counting and indexing loops
    {
      config.num_loops0 = 0;
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph);
	   vi != vi_end; ++vi) {
	if (config.bottom[*vi].loop_segment(0).index ==
	    loop_segment::undefined) {
	  if (config.bottom[*vi].loop_segment(0).root()->index ==
	      loop_segment::undefined)
	    config.bottom[*vi].loop_segment(0).root()->index =
	      (config.num_loops0)++;
	  config.bottom[*vi].loop_segment(0).index =
	    config.bottom[*vi].loop_segment(0).root()->index;
	}
	if (config.top[*vi].loop_segment(0).index ==
	    loop_segment::undefined) {
	  if (config.top[*vi].loop_segment(0).root()->index ==
	      loop_segment::undefined)
	    config.top[*vi].loop_segment(0).root()->index =
	      (config.num_loops0)++;
	  config.top[*vi].loop_segment(0).index =
	    config.top[*vi].loop_segment(0).root()->index;
	}
      }
    }
    {
      config.num_loops = config.num_loops0;
      operator_iterator oi_end = config.os.end();
      for (operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
	if (oi->loop_segment(0).index == loop_segment::undefined) {
	  if (oi->loop_segment(0).root()->index == loop_segment::undefined)
	    oi->loop_segment(0).root()->index = config.num_loops++;
	  oi->loop_segment(0).index = oi->loop_segment(0).root()->index;
	}
	if (oi->loop_segment(1).index == loop_segment::undefined) {
	  if (oi->loop_segment(1).root()->index == loop_segment::undefined)
	    oi->loop_segment(1).root()->index = config.num_loops++;
	  oi->loop_segment(1).index = oi->loop_segment(1).root()->index;
	}
      }
    }
  }
  template<class RNG>
  static void generate_loops(config_type& config, parameter_type& p,
			     RNG& uniform_01)
  {
    generate_loops(config, p.virtual_graph, p.model, p.beta, p.chooser,
		   uniform_01);
  }

  static double energy_offset(const vg_type& vg, const model_type& model)
  {
    double offset = 0.;
    typename alps::property_map<alps::bond_type_t, graph_type, int>::const_type
      bond_type(alps::get_or_default(alps::bond_type_t(), vg.graph, 0));
    edge_iterator ei_end = boost::edges(vg.graph).second;
    for (edge_iterator ei = boost::edges(vg.graph).first; ei != ei_end; ++ei)
      offset += model.bond(bond_type[*ei]).c() +
	weight_type(model.bond(bond_type[*ei])).offset();
    return offset;
  }
  static double energy_offset(const parameter_type& p)
  { return energy_offset(p.virtual_graph, p.model); }
}; // struct sse

} // end namespace looper

#endif // LOOPER_SSE_H
