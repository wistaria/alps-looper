/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

// $Id: sse.h 604 2004-01-16 08:35:21Z wistaria $

#ifndef LOOPER_SSE_H
#define LOOPER_SSE_H

#include <looper/permutation.h>
#include <looper/union_find.h>
#include <looper/virtual_graph.h>
#include <looper/weight.h>
#include <boost/throw_exception.hpp>
#include <cmath> // for std::sqrt
#include <stdexcept>

namespace looper {

typedef qmc_node sse_node;

template<class G, class M, class W = default_weight, class N = sse_node>
struct sse;

template<class G, class M, class W, class N>
struct sse<virtual_graph<G>, M, W, N>
{
  typedef virtual_graph<G>                      vg_type;
  typedef typename virtual_graph<G>::graph_type graph_type;
  typedef M                                     model_type;
  typedef W                                     weight_type;
  typedef N                                     node_type;

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
    typedef sse<virtual_graph<G>, M, W, N> qmc_type;

    typedef virtual_graph<G> vg_type;
    typedef typename vg_type::graph_type graph_type;
    typedef typename vg_type::mapping_type mapping_type;
    typedef M model_type;
    typedef W weight_type;

    template<class RG>
    parameter_type(const RG& rg, const model_type& m, double b)
      : virtual_graph(), model(m), beta(b), is_bipartite(false),
	chooser(), ez_offset(0.)
    { 
      looper::generate_virtual_graph(rg, model, virtual_graph);
      is_bipartite = alps::set_parity(virtual_graph.graph);
      chooser.init(virtual_graph, model);
      update_offset();
    }

    void update_offset()
    {
      typename boost::graph_traits<graph_type>::edge_iterator ei, ei_end;
      ez_offset = 0.;
      for (boost::tie(ei, ei_end) = boost::edges(virtual_graph.graph);
	   ei != ei_end; ++ei)
	ez_offset +=
	  weight_type(model.bond(bond_type(*ei, virtual_graph.graph))).
	    offset();
    }

    vg_type                   virtual_graph;
    const model_type&         model;
    double                    beta;
    bool                      is_bipartite;
    bond_chooser<weight_type> chooser;
    double                    ez_offset;
  };

  struct config_type
  {
    typedef sse<virtual_graph<G>, M, W, N> qmc_type;

    typedef N                                node_type;
    typedef std::vector<node_type>           os_type;
    typedef typename os_type::iterator       iterator;
    typedef typename os_type::const_iterator const_iterator;

    BOOST_STATIC_CONSTANT(bool, is_sse = true);
    BOOST_STATIC_CONSTANT(bool, is_path_integral = false);

    os_type bottom;
    os_type top;
    os_type os;
    unsigned int num_operators;

    unsigned int num_loops0;
    unsigned int num_loops;

    alps::ODump& save(alps::ODump& od) const {
      return od << bottom << top << os << num_operators;
    }
    alps::IDump& load(alps::IDump& id) {
      return id >> bottom >> top >> os >> num_operators;
    }
  };

  // helper functions: segment_d, segment_d

  static typename node_type::segment_type&
  segment_d(const typename config_type::iterator& itr, int leg)
  {
    if (itr->is_refl()) {
      return itr->loop_segment(1);
    } else {
      return itr->loop_segment(leg);
    }
  }

  static const typename node_type::segment_type&
  segment_d(const typename config_type::const_iterator& itr, int leg)
  {
    if (itr->is_refl()) {
      return itr->loop_segment(1);
    } else {
      return itr->loop_segment(leg);
    }
  }

  static typename node_type::segment_type&
  segment_u(const typename config_type::iterator& itr, int leg)
  {
    if (itr->is_refl()) {
      return itr->loop_segment(0);
    } else {
      return itr->loop_segment(1-leg);
    }
  }

  static const typename node_type::segment_type&
  segment_u(const typename config_type::const_iterator& itr, int leg)
  {
    if (itr->is_refl()) {
      return itr->loop_segment(0);
    } else {
      return itr->loop_segment(1-leg);
    }
  }

  //
  // update functions
  //
  
  // initialize
  static void initialize(config_type& config, const vg_type& vg, int ni = 16)
  {
    typedef typename config_type::iterator operator_iterator;

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

  static bool check_and_resize(config_type& config)
  {
    int old_size = config.os.size();
    if (config.num_operators > 0.75 * old_size) {
      config.os.resize(2 * old_size);
      for (int i = old_size - 1; i >= 0; --i) {
	config.os[2 * i + 1] = config.os[i];
	config.os[2 * i    ].clear_graph();
	config.os[2 * i    ].set_to_identity();
	config.os[    i    ].clear_graph();
	config.os[    i    ].set_to_identity();
      }
      return true;
    } else {
      return false;
    }
  }

  template<class RNG>
  static void generate_loops(config_type& config, const vg_type& vg,
			     const model_type& /* model */, double beta,
			     const bond_chooser<weight_type>& bc,
			     RNG& uniform_01)
  {
    typedef typename config_type::iterator operator_iterator;

    //
    // diagonal update & labeling
    //

    // check & resize
    check_and_resize(config);

    // copy spin configurations at the bottom
    std::vector<int> curr_conf(boost::num_vertices(vg.graph));
    {
      vertex_iterator vi, vi_end;
      std::vector<int>::iterator itr = curr_conf.begin();
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph); 
	   vi != vi_end; ++vi, ++itr) *itr = config.bottom[*vi].conf();
    }
    
    // scan over operators
    if (boost::num_edges(vg.graph) > 0) {
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
	if (!oi->is_identity()) {
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
  }
  template<class RNG>
  static void generate_loops(config_type& config, parameter_type& p,
			     RNG& uniform_01)
  {
    generate_loops(config, p.virtual_graph, p.model, p.beta, p.chooser,
		   uniform_01);
  }

  template<class RNG>
  static void flip_and_cleanup(config_type& config, const vg_type& vg,
                               RNG& uniform_01)
  {
    typedef typename config_type::iterator operator_iterator;
    typedef RNG rng_type;

    std::vector<int> flip(config.num_loops);
    std::generate(flip.begin(), flip.end(),
		  boost::variate_generator<rng_type&,
                                           boost::uniform_smallint<> >(
                    uniform_01, boost::uniform_smallint<>(0, 1)));

    // flip spins
    {
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph);
	   vi != vi_end; ++vi) {
	if (flip[config.bottom[*vi].loop_segment(0).index] == 1)
	  config.bottom[*vi].flip_conf();
	config.bottom[*vi].clear_graph();
	if (flip[config.top[*vi].loop_segment(0).index] == 1)
	  config.top[*vi].flip_conf();
	config.top[*vi].clear_graph();
      }
    }

    // update operators
    {
      operator_iterator oi_end = config.os.end();
      for (operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
	if ((!oi->is_identity()) &&
	    (flip[oi->loop_segment(0).index] ^
	     flip[oi->loop_segment(1).index] == 1))
	  oi->flip_operator();
	oi->clear_graph();
      }
    }
  }
  template<class RNG>
  static void flip_and_cleanup(config_type& config, const parameter_type& p,
                               RNG& uniform_01)
  {
    flip_and_cleanup(config, p.virtual_graph, uniform_01);
  }

  // measurements

  static int loop_index_0(int i, const config_type& config)
  { return config.bottom[i].loop_segment(0).index; }

  static int loop_index_1(int i, const config_type& config)
  { return config.top[i].loop_segment(0).index; }

  static double static_sz(int i, const config_type& config)
  { return 0.5 - (double)config.bottom[i].conf(); }

  static double dynamic_sz(int i, const config_type& config, const vg_type& vg)
  {
    typedef typename config_type::const_iterator const_operator_iterator;

    double sz = 0.;

    double c = static_sz(i, config);
    const_operator_iterator oi_end = config.os.end();
    for (const_operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
      sz += c;
      if (oi->is_offdiagonal()) {
	edge_iterator ei = boost::edges(vg.graph).first + oi->bond();
	if (boost::source(*ei, vg.graph) == i ||
	    boost::target(*ei, vg.graph) == i) c = -c;
      }
    }
    return sz / std::sqrt((double)config.os.size() *
    			  (double)(config.os.size() + 1));
  }

  static int num_offdiagonals(const config_type& config)
  {
    typedef typename config_type::const_iterator const_operator_iterator;
    int n = 0;
    const_operator_iterator oi_end = config.os.end();
    for (const_operator_iterator oi = config.os.begin(); oi != oi_end; ++oi)
      if (oi->is_offdiagonal()) ++n;
    return n;
  }

  static double energy_offset(const vg_type& vg, const model_type& model)
  {
    double offset = 0.;
    edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(vg.graph); ei != ei_end; ++ei)
      offset += model.bond(bond_type(*ei, vg.graph)).c();
    return offset;
  }

  // for debugging

  static void output(config_type& config, const parameter_type& p)
  {
    std::cout << "config at bottom: ";
    for (int i = 0; i < boost::num_vertices(p.graph); ++i)
      std::cout << (config.bottom[i].conf() == 0 ? '+' : '-');
    std::cout << std::endl;

    std::cout << "config at top   : ";
    for (int i = 0; i < boost::num_vertices(p.graph); ++i)
      std::cout << (config.top[i].conf() == 0 ? '+' : '-');
    std::cout << std::endl;

    std::cout << "operator string:\n";
    for (int i = 0; i < config.os.size(); ++i) {
      if (config.os[i].is_identity()) {
	std::cout << i << " identity\n";
      } else if (config.os[i].is_diagonal()) {
	std::cout << i << " diagonal     at bond " << config.os[i].bond()
		  << std::endl;
      } else {
	std::cout << i << " off-diagonal at bond " << config.os[i].bond()
		  << std::endl;
      }
    }
  }
  
}; // struct sse

} // end namespace looper

#endif // LOOPER_SSE_H
