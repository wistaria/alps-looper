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

// $Id: xxz.h 635 2004-02-29 10:59:51Z troyer $

#ifndef LOOPER_XXZ_H
#define LOOPER_XXZ_H

#include <alps/parameters.h>
#include <alps/lattice.h>
#include <alps/model.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace looper {

namespace detail {

template<class I>
I index(const alps::half_integer<I>& s0,
	const alps::half_integer<I>& s1,
	const alps::half_integer<I>& sz0,
	const alps::half_integer<I>& sz1)
{
  return s0.distance(sz0) * (s1.get_twice()+1) + s1.distance(sz1);
}

} // end namespace detail


class xxz_parameter
{
public:
  xxz_parameter() : c_(), jxy_(), jz_() {}
  xxz_parameter(double c, double jxy, double jz) : c_(c), jxy_(jxy), jz_(jz) {}

  double c() const { return c_; }
  double& c() { return c_; }
  
  double jxy() const { return jxy_; }
  double& jxy() { return jxy_; }

  double jz() const { return jz_; }
  double& jz() { return jz_; }

private:
  double c_;
  double jxy_;
  double jz_;
};


template <class T = double, class M = boost::numeric::ublas::matrix<T> >
class xxz_matrix
{
public:
  typedef T value_type;
  typedef M matrix_type;
  typedef typename matrix_type::size_type size_type;

  xxz_matrix() : matrix_() {}
  xxz_matrix(const xxz_matrix& m) : matrix_(m.matrix_) {}
  template<class I>
  xxz_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     double e0, double jxy, double jz) : matrix_()
  { build(s0, s1, e0, jxy, jz); }
  template<class I>
  xxz_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     const xxz_parameter& p) : matrix_()
  { build(s0, s1, p.c(), p.jxy(), p.jz()); }
  
  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  // access to to rows
  typename matrix_type::matrix_row_type operator[](size_type i) {
    return matrix_[i];
  }
  typename matrix_type::const_matrix_row_type operator[](size_type i) const {
    return matrix_[i];
  }

  template<class I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     double e0, double jxy, double jz)
  {
    typedef alps::half_integer<I> half_integer_type;
    
    // set matrix dimension
    int dim = (s0.get_twice()+1) * (s1.get_twice()+1);
    matrix_.resize(dim, dim);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
	matrix_[i][j] = value_type(0);

    // diagonal elements: jz sz0 sz1
    for (half_integer_type sz0 = s0; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1; --sz1) {
	matrix_[detail::index(s0, s1, sz0, sz1)]
	       [detail::index(s0, s1, sz0, sz1)] = 
	  e0 - jz * double(sz0) * double(sz1);
      }
    }

    // off-diagonal elements: jxy s0+ s1- / 2
    for (half_integer_type sz0 = s0-1; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1+1; --sz1) {
	matrix_[detail::index(s0, s1, sz0+1, sz1-1)]
               [detail::index(s0, s1, sz0, sz1)] =
	  - 0.5 * jxy *
	  std::sqrt(double(s0-sz0) * double(s0+sz0+1)) * 
	  std::sqrt(double(s1+sz1) * double(s1-sz1+1));
      }
    }

    // off-diagonal elements: jxy s0- s1+ / 2
    for (half_integer_type sz0 = s0; sz0 >= -s0+1; --sz0) {
      for (half_integer_type sz1 = s1-1; sz1 >= -s1; --sz1) {
	matrix_[detail::index(s0, s1, sz0-1, sz1+1)]
               [detail::index(s0, s1, sz0, sz1)] =
	  - 0.5 * jxy *
	  std::sqrt(double(s0+sz0) * double(s0-sz0+1)) * 
	  std::sqrt(double(s1-sz1) * double(s1+sz1+1));
      }
    }
  }

private:
  matrix_type matrix_;
};


template <class T, class M>
inline std::ostream& operator<<(std::ostream& os, const xxz_matrix<T, M>& m)
{
  os << m.matrix();
  return os;
}

//
// fitting a matrix to xxz_matrix
//

template <class I, class M>
inline boost::tuple<bool, typename M::value_type, typename M::value_type,
                    typename M::value_type>
fit2xxz(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1, 
	const M& mat, typename M::value_type tol = 1.0e-10)
{
  typedef M matrix_type;
  typedef typename M::value_type value_type;
  
  int dim = (s0.get_twice()+1) * (s1.get_twice()+1);
  xxz_matrix<value_type> m1(s0, s1, 0, 1, 1);
  
  // e0
  value_type e0 = 0.;
  for (int i = 0; i < dim; ++i) e0 += mat[i][i];
  e0 /= dim;
  
  // jz
  value_type jz = (mat[0][0] - e0) / m1[0][0];
  
  // jxy
  double jxy = 0;
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if ((i != j) && (m1[i][j] != 0)) {
	jxy = mat[i][j] / m1[i][j];
	break;
      }
    }
    if (jxy != 0) break;
  }
  
  // check
  bool success = true;
  xxz_matrix<value_type> m(s0, s1, e0, jxy, jz);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if (std::abs(mat[i][j] - m[i][j]) > tol) success = false;
    }
  }
  
  return boost::make_tuple(success, e0, jxy, jz);
}

template <class I, class T>
inline boost::tuple<bool, T, T, T>
fit2xxz(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1, 
	const boost::multi_array<T, 4>& mat, T tol = 1.0e-10)
{
  typedef T value_type;

  int d0 = s0.get_twice()+1;
  int d1 = s1.get_twice()+1;
  int dim = d0 * d1;

  boost::numeric::ublas::matrix<value_type> m(dim, dim);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
	for (int j1 = 0; j1 < d1; ++j1)
	  m[i0 * d1 + i1][j0 * d1 + j1] = mat[i0][i1][j0][j1];

  return fit2xxz(s0, s1, m, tol);
}


class xxz_model
{
public:
  typedef int type_type;
  typedef std::map<type_type, alps::half_integer<int> > spin_type;
  typedef std::map<type_type, xxz_parameter> bond_type;

  template<class G, class I>
  xxz_model(double Jxy, double Jz, const alps::half_integer<I>& spin,
	    const G& graph) : spin_(), bond_()
  { set_parameters(Jxy, Jz, spin, graph); }
  template<class G, class IntType>
  xxz_model(const alps::Parameters params, const G& graph,
	    const alps::ModelLibrary::OperatorDescriptorMap& ops,
	    const alps::HamiltonianDescriptor<IntType>& hd)
    : spin_(), bond_()
  { set_parameters(params, graph, ops, hd); }
  template<class G>
  xxz_model(const alps::Parameters params, const G& graph,
	    const alps::ModelLibrary& models) : spin_(), bond_()
  { set_parameters(params, graph, models); }
  
  template<class G, class I>
  void set_parameters(double Jxy, double Jz, const alps::half_integer<I>& spin,
		      const G& graph)
  {
    typedef G graph_type;
    typedef typename boost::graph_traits<graph_type>::vertex_iterator
      vertex_iterator;
    typedef typename boost::graph_traits<graph_type>::edge_iterator
      edge_iterator;
    
    // set site parameters
    typename alps::property_map<alps::site_type_t, graph_type,
                                type_type>::const_type
      site_type(alps::get_or_default(alps::site_type_t(), graph, 0));

    vertex_iterator vi_end = boost::vertices(graph).second;
    for (vertex_iterator vi = boost::vertices(graph).first; vi != vi_end;
	 ++vi) {
      type_type t = site_type[*vi];
      if (!spin_.count(t)) spin_[t] = spin;
    }

    // set bond parameters
    typename alps::property_map<alps::bond_type_t, graph_type,
                                type_type>::const_type
      bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));

    edge_iterator ei_end = boost::edges(graph).second;
    for (edge_iterator ei = boost::edges(graph).first; ei != ei_end; ++ei) {
      type_type t = bond_type[*ei];
      if (!bond_.count(t)) bond_[t] = xxz_parameter(0., Jxy, Jz);
    }
  }
  
  template<class G, class IntType>
  void set_parameters(alps::Parameters params, const G& graph,
		      const alps::ModelLibrary::OperatorDescriptorMap& ops,
		      const alps::HamiltonianDescriptor<IntType>& hd)
  {
    typedef G graph_type;
    typedef typename boost::graph_traits<graph_type>::vertex_iterator
      vertex_iterator;
    typedef typename boost::graph_traits<graph_type>::edge_iterator
      edge_iterator;
    
	// get default couplings
	params.copy_undefined(hd.default_parameters());
	
    // get site parameters
    typename alps::property_map<alps::site_type_t, graph_type,
                                type_type>::const_type
      site_type(alps::get_or_default(alps::site_type_t(), graph, 0));

    vertex_iterator vi_end = boost::vertices(graph).second;
    for (vertex_iterator vi = boost::vertices(graph).first; vi != vi_end;
	 ++vi) {
      type_type t = site_type[*vi];
      if (!spin_.count(t))
	spin_[t] = (double(hd.basis().site_basis(t).num_states()) - 1) / 2;
    }

    // get bond parameters
    std::map<boost::tuple<type_type, type_type, type_type>, bool> bond_visited;
    typename alps::property_map<alps::bond_type_t, graph_type,
                                type_type>::const_type
      bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));

    edge_iterator ei_end = boost::edges(graph).second;
    for (edge_iterator ei = boost::edges(graph).first; ei != ei_end; ++ei) {
      type_type bt = bond_type[*ei];
      type_type st0 = site_type[boost::source(*ei, graph)];
      type_type st1 = site_type[boost::target(*ei, graph)];
      if (!bond_visited[boost::make_tuple(bt, st0, st1)]) {
	bond_visited[boost::make_tuple(bt, st0, st1)] = true;
	boost::multi_array<double ,4> bm =
	  hd.bond_term(bt).
            template matrix<double>(hd.basis().site_basis(st0),
				    hd.basis().site_basis(st1),
				    ops, params);
	boost::tuple<bool, double, double, double>
	  fit = fit2xxz(spin_[st0], spin_[st1], bm);
	if (!fit.template get<0>())
	  boost::throw_exception(std::runtime_error("fitting to XXZ model failed"));
	if (!bond_.count(bt)) {
	  bond_[bt] = xxz_parameter(fit.template get<1>(),
				    fit.template get<2>(),
				    fit.template get<3>());
	} else {
	  if (bond_[bt].c() != fit.template get<1>() ||
	      bond_[bt].jxy() != fit.template get<2>() ||
	      bond_[bt].jz() != fit.template get<3>())
	    boost::throw_exception(std::runtime_error("inconsistent bond parameter(s)"));
	}
      }
    }
  }
  
  template<class G>
  void set_parameters(const alps::Parameters params,
		      const G& graph,
		      const alps::ModelLibrary& models)
  {
    // get Hamilton operator from ModelLibrary
    alps::HamiltonianDescriptor<short> hd(models.hamiltonian(params["MODEL"]));
    alps::Parameters p(params);
    p.copy_undefined(hd.default_parameters());
    hd.set_parameters(p);
    set_parameters(p, graph, models.simple_operators(), hd);
  }
  
  int num_spin_types() const { return spin_.size(); }
  bool is_uniform_spin() const { return num_spin_types() == 1; }
  alps::half_integer<int> spin(type_type t) const
    { return spin_.find(t)->second; }
  alps::half_integer<int> uniform_spin() const
  {
    assert(is_uniform_spin());
    return spin_.begin()->second;
  }

  int num_bond_types() const { return bond_.size(); }
  bool is_uniform_bond() const { return num_bond_types() == 1; }
  xxz_parameter bond(type_type t) const
    { return bond_.find(t)->second; }
  xxz_parameter uniform_bond() const
  {
    assert(is_uniform_bond());
    return bond_.begin()->second;
  }

private:
  std::map<type_type, alps::half_integer<int> > spin_;
  std::map<type_type, xxz_parameter> bond_;
};

} // end namespace looper

#endif // LOOPER_XXZ_H
