/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_XXZ_H
#define LOOPER_XXZ_H

#include <looper/lapack.h>
#include <looper/operator.h>
#include <looper/util.h>

#include <alps/parameters.h>
#include <alps/lattice.h>
#include <alps/model.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace looper {

//
// forward declarations
//

template<typename I = int> class site_parameter;

class xxz_parameter;

template<typename T = double> class site_matrix;

template<typename T = double> class xxz_matrix;

template<typename I = int> class xxz_model;

template<typename T, typename I> bool fit2site(const boost::multi_array<T, 2>& mat, site_parameter<I>& param, T tol = 1.0e-10);

template<typename T, typename I> bool fit2site(const site_matrix<T>& mat, site_parameter<I>& param, T tol = 1.0e-10);

template<typename T> bool fit2xxz(const boost::multi_array<T, 4>& mat, xxz_parameter& param, T tol = 1.0e-10);

template<typename T> bool fit2xxz(const xxz_matrix<T>& mat, xxz_parameter& param, T tol = 1.0e-10);


//
// paramters
//

template<typename I>
class site_parameter
{
public:
  typedef alps::half_integer<I> spin_type;

  site_parameter() : s_(0.5), c_(0), hx_(0), hz_(0) {}
  site_parameter(double s) : s_(s), c_(0), hx_(0), hz_(0) {}
  template<typename J>
  site_parameter(const alps::half_integer<J>& s) :
    s_(s), c_(0), hx_(0), hz_(0) {}
  site_parameter(double s, double c, double hx, double hz) :
    s_(s), c_(c), hx_(hx), hz_(hz) {}
  template<typename J>
  site_parameter(const alps::half_integer<J>& s,
		 double c, double hx, double hz) :
    s_(s), c_(c), hx_(hx), hz_(hz) {}
  site_parameter(const boost::multi_array<double, 2>& mat)
  {
    bool success = fit2site(mat, *this);
    if (!success)
      boost::throw_exception(std::runtime_error("fitting to site_parameter failed"));
  }

  template<typename J>
  bool operator==(const site_parameter<J>& rhs) const
  {
    return nearly_equal(s(), rhs.s()) && nearly_equal(c(), rhs.c()) &&
      nearly_equal(hx(), rhs.hx()) && nearly_equal(hz(), rhs.hz());
  }

  template<typename J>
  bool operator!=(const site_parameter<J>& rhs) const
  {
    return !operator==(rhs);
  }

  const spin_type& s() const { return s_; }
  spin_type& s() { return s_; }

  double c() const { return c_; }
  double& c() { return c_; }

  double hx() const { return hx_; }
  double& hx() { return hx_; }

  double hz() const { return hz_; }
  double& hz() { return hz_; }


private:
  spin_type s_;
  double c_, hx_, hz_;
};


class xxz_parameter
{
public:
  xxz_parameter() : c_(), jxy_(), jz_() {}
  xxz_parameter(double c, double jxy, double jz) : c_(c), jxy_(jxy), jz_(jz) {}
  xxz_parameter(const boost::multi_array<double, 4>& mat)
  {
    bool success = fit2xxz(mat, *this);
    if (!success)
      boost::throw_exception(std::runtime_error("fitting to xxz_parameter failed"));
  }

  double c() const { return c_; }
  double& c() { return c_; }

  double jxy() const { return jxy_; }
  double& jxy() { return jxy_; }

  double jz() const { return jz_; }
  double& jz() { return jz_; }

  bool operator==(const xxz_parameter& rhs) const
  {
    return nearly_equal(c(), rhs.c()) && nearly_equal(jxy(), rhs.jxy()) &&
      nearly_equal(jz(), rhs.jz());
  }
  bool operator!=(const xxz_parameter& rhs) const
  {
    return !operator==(rhs);
  }

private:
  double c_;
  double jxy_;
  double jz_;
};


//
// matrices
//

template<typename T>
class site_matrix
{
public:
  typedef T value_type;
  typedef boost::multi_array<value_type, 2> matrix_type;
  typedef typename matrix_type::size_type size_type;

  site_matrix() : matrix_() {}
  site_matrix(const site_matrix& m) : matrix_(m.matrix_) {}
  template<typename I>
  site_matrix(const alps::half_integer<I>& s, value_type c, value_type hx,
	      value_type hz) : matrix_()
  { build(s, c, hx, hz); }
  template<typename I>
  site_matrix(const site_parameter<I>& s) : matrix_()
  { build(s.s(), s.c(), s.hx(), s.hz()); }

  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  template<typename I>
  void build(const alps::half_integer<I>& s, value_type c, value_type hx,
	     value_type hz)
  {
    typedef I                                integer_type;
    typedef alps::half_integer<integer_type> half_integer_type;

    // set matrix dimension
    int dim = s.get_twice()+1;
    matrix_.resize(boost::extents[dim][dim]);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
	matrix_[i][j] = value_type(0);

    // diagonal elements: c - hz sz
    for (half_integer_type sz = -s; sz <= s; ++sz)
      matrix_[sz.distance(-s)][sz.distance(-s)] = c - hz * value_type(sz);

    // off-diagonal elements: hx s+ / 2
    for (half_integer_type sz = -s; sz <= s-1; ++sz)
      matrix_[sz.distance(-s)+1][sz.distance(-s)] =
	- 0.5 * hx * std::sqrt(value_type(s-sz) * value_type(s+sz+1));

    // off-diagonal elements: hx s- / 2
    for (half_integer_type sz = -s+1; sz <= s; ++sz)
      matrix_[sz.distance(-s)-1][sz.distance(-s)] = 
	- 0.5 * hx * std::sqrt(value_type(s+sz) * value_type(s-sz+1));
  }

private:
  matrix_type matrix_;
};


template<typename T>
class xxz_matrix
{
public:
  typedef T value_type;
  typedef boost::multi_array<value_type, 4> matrix_type;
  typedef typename matrix_type::size_type size_type;

  xxz_matrix() : matrix_() {}
  xxz_matrix(const xxz_matrix& m) : matrix_(m.matrix_) {}
  template<typename I>
  xxz_matrix(const site_parameter<I>& s0, const site_parameter<I>& s1,
             value_type e0, value_type jxy, value_type jz) : matrix_()
  { build(s0.s(), s1.s(), e0, jxy, jz); }
  template<typename I>
  xxz_matrix(const site_parameter<I>& s0, const site_parameter<I>& s1,
             const xxz_parameter& p) : matrix_()
  { build(s0.s(), s1.s(), p.c(), p.jxy(), p.jz()); }

  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  template<typename I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
             value_type c, value_type jxy, value_type jz)
  {
    typedef I                                integer_type;
    typedef alps::half_integer<integer_type> half_integer_type;

    // set matrix dimension
    int d0 = s0.get_twice()+1;
    int d1 = s1.get_twice()+1;
    matrix_.resize(boost::extents[d0][d1][d0][d1]);
    for (int i0 = 0; i0 < d0; ++i0)
      for (int i1 = 0; i1 < d1; ++i1)
	for (int j0 = 0; j0 < d0; ++j0)
	  for (int j1 = 0; j1 < d1; ++j1)
	    matrix_[i0][i1][j0][j1] = value_type(0);

    // diagonal elements: c - jz sz0 sz1
    for (half_integer_type sz0 = -s0; sz0 <= s0; ++sz0) {
      for (half_integer_type sz1 = -s1; sz1 <= s1; ++sz1) {
	matrix_[sz0.distance(-s0)][sz1.distance(-s1)]
	  [sz0.distance(-s0)][sz1.distance(-s1)] =
	  c - jz * value_type(sz0) * value_type(sz1);
      }
    }

    // off-diagonal elements: - jxy s0+ s1- / 2
    for (half_integer_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
      for (half_integer_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
	matrix_[sz0.distance(-s0)+1][sz1.distance(-s1)-1]
	  [sz0.distance(-s0)][sz1.distance(-s1)] =
	  - 0.5 * jxy *
	  std::sqrt(value_type(s0-sz0) * value_type(s0+sz0+1)) *
	  std::sqrt(value_type(s1+sz1) * value_type(s1-sz1+1));
      }
    }

    // off-diagonal elements: - jxy s0- s1+ / 2
    for (half_integer_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
      for (half_integer_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
	matrix_[sz0.distance(-s0)-1][sz1.distance(-s1)+1]
	  [sz0.distance(-s0)][sz1.distance(-s1)] =
	  - 0.5 * jxy *
	  std::sqrt(value_type(s0+sz0) * value_type(s0-sz0+1)) *
	  std::sqrt(value_type(s1-sz1) * value_type(s1+sz1+1));
      }
    }
  }

private:
  matrix_type matrix_;
};


//
// models
//

template<typename I>
class xxz_model
{
public:
  typedef int type_type;
  typedef std::map<type_type, site_parameter<I> > site_type;
  typedef std::map<type_type, xxz_parameter> bond_type;
  typedef typename site_parameter<I>::spin_type spin_type;

  template<typename IntType, typename G>
  xxz_model(double Jxy, double Jz, const alps::half_integer<IntType>& spin,
            const G& graph) : site_(), bond_()
  { set_parameters(Jxy, Jz, spin, graph); }
  template<typename G, typename IntType>
  xxz_model(const alps::Parameters params, const G& graph,
            const alps::ModelLibrary::OperatorDescriptorMap& ops,
            const alps::HamiltonianDescriptor<IntType>& hd)
    : site_(), bond_()
  { set_parameters(params, graph, ops, hd); }
  template<typename G>
  xxz_model(const alps::Parameters params, const G& graph,
            const alps::ModelLibrary& models) : site_(), bond_()
  { set_parameters(params, graph, models); }

  template<typename IntType, typename G>
  void set_parameters(double Jxy, double Jz,
                      const alps::half_integer<IntType>& spin,
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
      if (!site_.count(t)) site_[t].s() = spin;
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

  template<typename G, typename IntType>
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
      if (!site_.count(t))
        site_[t] = site_parameter<I>(alps::get_matrix(0.,hd.site_term(t),
          hd.basis().site_basis(t), ops, params));
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
        xxz_parameter p(alps::get_matrix(0., hd.bond_term(bt), 
					 hd.basis().site_basis(st0),
					 hd.basis().site_basis(st1),
					 ops, params));
        if (!bond_.count(bt)) {
          bond_[bt] = p;
        } else {
          if (bond_[bt] != p)
            boost::throw_exception(std::runtime_error("inconsistent bond parameter(s)"));
        }
      }
    }
  }

  template<typename G>
  void set_parameters(const alps::Parameters params,
                      const G& graph,
                      const alps::ModelLibrary& models)
  {
    // get Hamilton operator from ModelLibrary
    alps::HamiltonianDescriptor<short>
      hd(models.get_hamiltonian(params["MODEL"]));
    alps::Parameters p(params);
    p.copy_undefined(hd.default_parameters());
    hd.set_parameters(p);
    set_parameters(p, graph, models.operators(), hd);
  }

  int num_site_types() const { return site_.size(); }
  bool is_uniform_site() const { return num_site_types() == 1; }
  site_parameter<I> site(type_type t) const
  { return site_.find(t)->second; }
  site_parameter<I> uniform_site() const
  {
#ifndef NDEBUG
    assert(is_uniform_site());
#endif
    return site_.begin()->second;
  }

  int num_bond_types() const { return bond_.size(); }
  bool is_uniform_bond() const { return num_bond_types() == 1; }
  xxz_parameter bond(type_type t) const
    { return bond_.find(t)->second; }
  xxz_parameter uniform_bond() const
  {
#ifndef NDEBUG
    assert(is_uniform_bond());
#endif
    return bond_.begin()->second;
  }

private:
  std::map<type_type, site_parameter<I> > site_;
  std::map<type_type, xxz_parameter> bond_;
};


//
// functions for fitting
//

template<typename T, typename I>
bool fit2site(const boost::multi_array<T, 2>& mat, site_parameter<I>& param, T tol)
{
  typedef T value_type;

#ifndef NDEBUG
  assert(mat.shape()[0] == mat.shape()[1]);
#endif

  int dim = mat.shape()[0];
  int m = dim * dim;
  int n = 3;
  
  typename site_parameter<I>::spin_type s((double)(dim - 1)/2);

  alps::SiteBasisDescriptor<short> basis(spin_basis(s));
  boost::multi_array<value_type, 2> mat_c(
    alps::get_matrix(value_type(),
		     alps::SiteTermDescriptor<short>("1"),
		     basis, spin_operators(), alps::Parameters()));
  boost::multi_array<value_type, 2> mat_hx(
    alps::get_matrix(value_type(),
		     alps::SiteTermDescriptor<short>("-(Splus+Sminus)/2"),
		     basis, spin_operators(), alps::Parameters()));
  boost::multi_array<value_type, 2> mat_hz(
    alps::get_matrix(value_type(),
		     alps::SiteTermDescriptor<short>("-Sz"),
		     basis, spin_operators(), alps::Parameters()));

  boost::numeric::ublas::matrix<value_type,
    boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<value_type> b(m);
  boost::numeric::ublas::vector<value_type> x(n);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      int k = dim * i + j;
      b(k) = mat[i][j];
      a(k, 0) = mat_c[i][j];
      a(k, 1) = mat_hx[i][j];
      a(k, 2) = mat_hz[i][j];
    }
  }
  
  // call linear least sqaure problem solver
  bool success = (solve_llsp(a, b, x) < tol);

  value_type c = numeric_cast<value_type>(x(0));
  if (std::abs(c) < tol) c = 0;
  value_type hx = numeric_cast<value_type>(x(1));
  if (std::abs(hx) < tol) hx = 0;
  value_type hz = numeric_cast<value_type>(x(2));
  if (std::abs(hz) < tol) hz = 0;

  if (success) param = site_parameter<I>(s, c, hx, hz);

  return success;
}


template<typename T, typename I>
bool fit2site(const site_matrix<T>& mat, site_parameter<I>& param, T tol)
{
  return fit2site(mat.matrix(), param, tol);
}


template<typename T>
bool fit2xxz(const boost::multi_array<T, 4>& mat, xxz_parameter& param, T tol)
{
#ifndef NDEBUG
  assert(mat.shape()[0] == mat.shape()[2]);
  assert(mat.shape()[1] == mat.shape()[3]);
#endif

  typedef T value_type;

  int d0 = mat.shape()[0];
  int d1 = mat.shape()[1];
  int dim = d0 * d1;
  int m = dim * dim;
  int n = 3;

  alps::Parameters p;
  p["S0"] = (double)(d0 - 1)/2;
  p["S1"] = (double)(d1 - 1)/2;
  alps::SiteBasisDescriptor<short> basis0(spin_basis((double)(d0 - 1)/2));
  alps::SiteBasisDescriptor<short> basis1(spin_basis((double)(d1 - 1)/2));

  std::set<std::string> suffixes;
  suffixes.insert("0");
  suffixes.insert("1");

  boost::multi_array<value_type, 4> mat_c(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor<short>("1"),
                     basis0, basis1, spin_operators(suffixes),
		     p));
  boost::multi_array<value_type, 4> mat_jxy(
    alps::get_matrix(value_type(),
		     alps::BondTermDescriptor<short>(
		       "-(Splus0(i)*Sminus1(j)+Sminus0(i)*Splus1(j))/2"),
                     basis0, basis1, spin_operators(suffixes),
		     p));
  boost::multi_array<value_type, 4> mat_jz(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor<short>("-Sz0(i)*Sz1(j)"),
                     basis0, basis1, spin_operators(suffixes), 
		     p));

  boost::numeric::ublas::matrix<value_type,
    boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<value_type> b(m);
  boost::numeric::ublas::vector<value_type> x(n);
  for (int i0 = 0; i0 < d0; ++i0) {
    for (int i1 = 0; i1 < d1; ++i1) {
      for (int j0 = 0; j0 < d0; ++j0) {
        for (int j1 = 0; j1 < d1; ++j1) {
          int k = dim * (i0 * d1 + i1) + (j0 * d1 + j1);
          b(k) = mat[i0][i1][j0][j1];
          a(k, 0) = mat_c[i0][i1][j0][j1];
          a(k, 1) = mat_jxy[i0][i1][j0][j1];
          a(k, 2) = mat_jz[i0][i1][j0][j1];
        }
      }
    }
  }

  // call linear least sqaure problem solver
  bool success = (solve_llsp(a, b, x) < tol);

  value_type c = numeric_cast<value_type>(x(0));
  if (std::abs(c) < tol) c = 0;
  value_type jxy = numeric_cast<value_type>(x(1));
  if (std::abs(jxy) < tol) jxy = 0;
  value_type jz = numeric_cast<value_type>(x(2));
  if (std::abs(jz) < tol) jz = 0;

  if (success) param = xxz_parameter(c, jxy, jz);

  return success;
}


template<typename T>
bool fit2xxz(const xxz_matrix<T>& mat, xxz_parameter& param, T tol)
{
  return fit2xxz(mat.matrix(), param, tol);
}


} // end namespace looper


#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<typename T>
std::ostream& operator<<(std::ostream& os, const looper::site_matrix<T>& m)
{
  os << m.matrix();
  return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const looper::xxz_matrix<T>& m)
{
  boost::numeric::ublas::matrix<T> mat;
  flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

template<typename I>
std::ostream& operator<<(std::ostream& os, const looper::site_parameter<I>& p)
{
  os << "C = " << p.c() << ", Hx = " << p.hx() << ", Hz = " << p.hz();
  return os;
}

std::ostream& operator<<(std::ostream& os, const looper::xxz_parameter& p)
{
  os << "C = " << p.c() << ", Jxy = " << p.jxy() << ", Jz = " << p.jz();
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace looper
#endif

#endif // LOOPER_XXZ_H
