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

// $Id: exact_diag.h 561 2003-11-12 15:53:43Z wistaria $

#ifndef LOOPER_EXACT_DIAG_H
#define LOOPER_EXACT_DIAG_H

#include <looper/lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace looper {

template<class G, class M, class T = double>
struct exact_diagonalization
{
  typedef G graph_type;
  typedef M model_type;
  typedef T value_type;
  typedef boost::numeric::ublas::vector<value_type> vector_type;
  typedef boost::numeric::ublas::matrix<value_type> matrix_type;

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
    parameter_type(graph_type& g, const model_type& m, double b)
      : graph(g), model(m), beta(b)
    {
      is_bipartite = alps::set_parity(graph);
    }

    graph_type& graph;
    const model_type& model;
    double            beta;
    bool              is_bipartite;
  };

  struct config_type
  {
    int dimension;
    std::vector<std::pair<int, int> > basis;
    matrix_type hamiltonian;
    vector_type eigenvalues;
  };

  // helper functions

  static int index(int s0, int s1, int /* d0 */, int d1)
  { return s0 * d1 + s1; }
  static int sz2(int s, int i, const config_type& config)
  { return (s / config.basis[i].second) % config.basis[i].first; }
  static int up(int s, int i, const config_type& config)
  {
    return (s > 0 && sz2(s, i, config) != config.basis[i].first - 1) ?
      (s + config.basis[i].second) : -1;
  }
  static int down(int s, int i, const config_type& config)
  {
    return (s > 0 && sz2(s, i, config) != 0) ?
      (s - config.basis[i].second) : -1;
  }

  static generate_matrix(const parameter_type& param, config_type& config)
  {
    // setup basis
    config.dimension = 0;
    config.basis.clear();
    vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.graph);
	 vi != vi_end; ++vi) {
      int d = param.model.spin(site_type(*vi, param.graph)).get_twice() + 1;
      if (config.basis.size() == 0) config.dimension = 1;
      config.basis.push_back(std::make_pair(d, config.dimension));
      config.dimension *= d;
    }

    // fill Hamiltonian matrix
    config.hamiltonian.clear();
    config.hamiltonian.resize(config.dimension, config.dimension);
    edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(param.graph);
	 ei != ei_end; ++ei) {
      vertex_descriptor v0 = boost::source(*ei, param.graph);
      vertex_descriptor v1 = boost::target(*ei, param.graph);
      int d0 = config.basis[v0].first;
      int d1 = config.basis[v1].first;

      xxz_matrix<value_type, matrix_type>
	m(param.model.spin(site_type(v0, param.graph)),
	  param.model.spin(site_type(v1, param.graph)),
	  param.model.bond(bond_type(*ei, param.graph)));

      for (int s = 0; s < config.dimension; ++s) {
	int s0 = sz2(s, v0, config);
	int s1 = sz2(s, v1, config);
	{
	  // diagonal element
	  config.hamiltonian(s, s) +=
	    m[index(s0, s1, d0, d1)][index(s0, s1, d0, d1)];
	}
	{
	  // up-down
	  int t = down(up(s, v0, config), v1, config);
	  if (t > 0) {
	    int t0 = sz2(t, v0, config);
	    int t1 = sz2(t, v1, config);
	    config.hamiltonian(s, t) +=
	      m[index(s0, s1, d0, d1)][index(t0, t1, d0, d1)];
	  }
	}
	{
	  // down-up
	  int t = up(down(s, v0, config), v1, config);
	  if (t > 0) {
	    int t0 = sz2(t, v0, config);
	    int t1 = sz2(t, v1, config);
	    config.hamiltonian(s, t) +=
	      m[index(s0, s1, d0, d1)][index(t0, t1, d0, d1)];
	  }
	}
      }
    }

    config.eigenvalues.resize(config.dimension);
  }

  static void diagonalize(const parameter_type& param, config_type& config)
  {
    looper::diagonalize(config.hamiltonian, config.eigenvalues, true);
  }

  static boost::tuple<double, double, double>
  energy(const parameter_type& param, const config_type& config)
  {
    double part = 0.;
    double ene = 0.;
    double ene2 = 0.;
    double gs_ene = config.eigenvalues(0);
    std::cout << gs_ene << std::endl;
    typename vector_type::const_reverse_iterator ev_end =
      config.eigenvalues.rend();
    for (typename vector_type::const_reverse_iterator ev =
	   config.eigenvalues.rbegin(); ev != ev_end; ++ev) {
      std::cout << *ev << std::endl;
      // Boltzman weight
      double weight;
      weight = std::exp(-param.beta * (*ev - gs_ene));

      // partition function
      part += weight;
    
      // energy
      ene += (*ev) * weight;
      ene2 += (*ev) * (*ev) * weight;
    }
    ene /= part;
    ene2 /= part;
    double c = param.beta * param.beta * (ene2 - ene * ene);

    ene /= boost::num_vertices(param.graph);
    ene2 /= boost::num_vertices(param.graph);
    c /= boost::num_vertices(param.graph);
    
    boost::make_tuple(ene, ene2, c);
  }
};

} // end namespace looper

#endif // LOOPER_EXACT_DIAG_H
