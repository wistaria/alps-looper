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

// $Id: exact_diag.h 564 2003-11-13 08:08:07Z wistaria $

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
    template<class RG>
    parameter_type(const RG& g, const model_type& m, double b)
      : graph(), model(m), beta(b)
    {
      alps::copy_graph(g, graph);
      is_bipartite = alps::set_parity(graph);
    }

    graph_type        graph;
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
    mutable vector_type mtab;
    mutable vector_type stab;
  };

  // helper functions

  static int index(int s0, int s1, int /* d0 */, int d1)
  { return s0 * d1 + s1; }
  static int sz2(int s, int i, const config_type& config)
  { return (s / config.basis[i].second) % config.basis[i].first; }
  static double sz(int s, int i, const config_type& config)
  { return ((double)config.basis[i].first - 1) / 2 - sz2(s, i, config); }
  static int up(int s, int i, const config_type& config)
  {
    return (s >= 0 && sz2(s, i, config) != config.basis[i].first - 1) ?
      (s + config.basis[i].second) : -1;
  }
  static int down(int s, int i, const config_type& config)
  {
    return (s >= 0 && sz2(s, i, config) != 0) ?
      (s - config.basis[i].second) : -1;
  }

  static void generate_mtab(const parameter_type& param,
			    const config_type& config)
  {
    if (config.mtab.size() == 0) {
      config.mtab.resize(config.dimension);
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(param.graph);
	   vi != vi_end; ++vi) {
	typename vector_type::iterator s = config.mtab.begin();
	typename vector_type::iterator s_end = config.mtab.end();
	int state = 0;
	for (; s != s_end; ++s, ++state) {
	  *s += sz(state, *vi, config);
	}
      }
    }
  }

  static void generate_stab(const parameter_type& param, 
			    const config_type& config)
  {
    if (param.is_bipartite && config.stab.size() == 0) {
      config.stab.resize(config.dimension);
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(param.graph);
	   vi != vi_end; ++vi) {
	typename vector_type::iterator s = config.stab.begin();
	typename vector_type::iterator s_end = config.stab.end();
	int state = 0;
	for (; s != s_end; ++s, ++state) {
	  *s += gauge(*vi, param.graph) * sz(state, *vi, config);
	}
      }
    }
  }

  static void generate_matrix(const parameter_type& param, config_type& config)
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

  static void diagonalize(const parameter_type& /* param */,
			  config_type& config)
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
    typename vector_type::const_reverse_iterator ev_end =
      config.eigenvalues.rend();
    for (typename vector_type::const_reverse_iterator ev =
	   config.eigenvalues.rbegin(); ev != ev_end; ++ev) {
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
    
    return boost::make_tuple(ene, ene2, c);
  }

  static boost::tuple<double, double, double, double>
  magnetization(const parameter_type& param, const config_type& config)
  {
    generate_mtab(param, config);
    generate_stab(param, config);

    double part = 0.;
    double gs_ene = config.eigenvalues(0);
    double umag2 = 0.;
    double smag2 = 0.;
    double usus = 0.;
    double ssus = 0.;

    typename vector_type::const_reverse_iterator ev =
      config.eigenvalues.rbegin();
    typename vector_type::const_reverse_iterator ev_end =
      config.eigenvalues.rend();
    typename matrix_type::const_reverse_iterator1 evec =
      config.hamiltonian.rbegin1();
    for (; ev != ev_end; ++ev, ++evec) {
      // Boltzman weight
      double weight;
      weight = std::exp(-param.beta * (*ev - gs_ene));

      // partition function
      part += weight;
    
      // uniform and staggered magnetization
      if (param.is_bipartite) {
	double m2 = 0.;
	double s2 = 0.;
	typename matrix_type::const_iterator2 j = evec.begin();
	typename vector_type::const_iterator itr_m = config.mtab.begin();
	typename vector_type::const_iterator itr_s = config.stab.begin();
	for (; j != evec.end(); ++j, ++itr_m, ++itr_s) {
	  double w = std::pow(*j, 2);
	  m2 += std::pow(*itr_m, 2) * w;
	  s2 += std::pow(*itr_s, 2) * w;
	}
	umag2 += m2 * weight;
	smag2 += s2 * weight;
      } else {
	double m2 = 0.;
	typename matrix_type::const_iterator2 j = evec.begin();
	typename vector_type::const_iterator itr_m = config.mtab.begin();
	typename vector_type::const_iterator itr_s = config.stab.begin();
	for (; j != evec.end(); ++j, ++itr_m, ++itr_s) {
	  double w = std::pow(*j, 2);
	  m2 += std::pow(*itr_m, 2) * w;
	}
	umag2 += m2 * weight;
      }

      // uniform and staggered susceptibility
      if (param.is_bipartite) {
	typename matrix_type::const_reverse_iterator1 evec_j =
	  config.hamiltonian.rbegin1();
	typename vector_type::const_reverse_iterator ev_j =
	  config.eigenvalues.rbegin();
	for (; evec_j != evec; ++evec_j, ++ev_j) {
	  // for evec_j != evec
	  double mu = 0.;
	  double su = 0.;
	  typename matrix_type::const_iterator2 j = evec_j.begin();
	  typename vector_type::const_iterator itr_m = config.mtab.begin();
	  typename vector_type::const_iterator itr_s = config.stab.begin();
	  for (typename matrix_type::const_iterator2 k = evec.begin();
	       k != evec.end(); ++k, ++j, ++itr_m, ++itr_s) {
	    mu += (*k) * (*itr_m) * (*j);
	    su += (*k) * (*itr_s) * (*j);
	  }
	  double wij;
	  if (std::abs(*ev - *ev_j) > 1.e-12) {
	    wij = - (weight - exp(-param.beta * (*ev_j - gs_ene)))
	      / (*ev - *ev_j);
	  } else {
	    wij = param.beta * weight;
	  }
	  usus += 2 * std::pow(mu, 2) * wij;
	  ssus += 2 * std::pow(su, 2) * wij;
	}
	{
	  // for evec_j = evec
	  double mu = 0.;
	  double su = 0.;
	  typename matrix_type::const_iterator2 j = evec_j.begin();
	  typename vector_type::const_iterator itr_m = config.mtab.begin();
	  typename vector_type::const_iterator itr_s = config.stab.begin();
	  for (typename matrix_type::const_iterator2 k = evec.begin();
	       k != evec.end(); ++k, ++j, ++itr_m, ++itr_s) {
	    mu += (*k) * (*itr_m) * (*j);
	    su += (*k) * (*itr_s) * (*j);
	  }
	  double wij = param.beta * weight;
	  usus += std::pow(mu, 2) * wij;
	  ssus += std::pow(su, 2) * wij;
	}
      } else {
	typename matrix_type::const_reverse_iterator1 evec_j =
	  config.hamiltonian.rbegin1();
	typename vector_type::const_reverse_iterator ev_j =
	  config.eigenvalues.rbegin();
	for (; evec_j != evec; ++evec_j, ++ev_j) {
	  // for evec_j != evec
	  double mu = 0.;
	  typename matrix_type::const_iterator2 j = evec_j.begin();
	  typename vector_type::const_iterator itr_m = config.mtab.begin();
	  for (typename matrix_type::const_iterator2 k = evec.begin();
	       k != evec.end(); ++k, ++j, ++itr_m) {
	    mu += (*k) * (*itr_m) * (*j);
	  }
	  double wij;
	  if (std::abs(*ev - *ev_j) > 1.e-12) {
	    wij = - (weight - exp(-param.beta * (*ev_j - gs_ene)))
	      / (*ev - *ev_j);
	  } else {
	    wij = param.beta * weight;
	  }
	  usus += 2 * std::pow(mu, 2) * wij;
	}
	{
	  // for evec_j = evec
	  double mu = 0.;
	  typename matrix_type::const_iterator2 j = evec_j.begin();
	  typename vector_type::const_iterator itr_m = config.mtab.begin();
	  for (typename matrix_type::const_iterator2 k = evec.begin();
	       k != evec.end(); ++k, ++j, ++itr_m) {
	    mu += (*k) * (*itr_m) * (*j);
	  }
	  double wij = param.beta * weight;
	  usus += std::pow(mu, 2) * wij;
	}
      }
    }

    umag2 /= part;
    usus /= part;
    smag2 /= part;
    ssus /= part;

    umag2 /= boost::num_vertices(param.graph);
    usus /= boost::num_vertices(param.graph);
    smag2 /= boost::num_vertices(param.graph);
    ssus /= boost::num_vertices(param.graph);

    return boost::make_tuple(umag2, usus, smag2, ssus);
  }
};

} // end namespace looper

#endif // LOOPER_EXACT_DIAG_H
