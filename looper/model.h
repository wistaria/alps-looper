/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_MODEL_H
#define LOOPER_MODEL_H

#include <looper/graph.h>
#include <looper/lapack.h>
#include <looper/operator.h>
#include <looper/util.h>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace looper {

//
// forward declarations
//

// site_parameter
class site_parameter_hxz; // most generic
class site_parameter_noh; // default

// bond_parameter
class bond_parameter_xxz; // most generic & default

template<typename T = double> class site_matrix; // support hxz
template<typename T = double> class bond_matrix; // support xxz

template<typename T> bool fit2site(const boost::multi_array<T, 2>& mat,
  site_parameter_hxz& param, T tol = 1.0e-10);
template<typename T> bool fit2bond(const boost::multi_array<T, 4>& mat,
  bond_parameter_xxz& param, T tol = 1.0e-10);

template<typename SITE_P = site_parameter_noh,
         typename BOND_P = bond_parameter_xxz>
class model_parameter;


//
// parameters
//

class site_parameter_hxz
{
public:
  typedef alps::half_integer<int> spin_type;

  BOOST_STATIC_CONSTANT(bool, has_hx = true);
  BOOST_STATIC_CONSTANT(bool, has_hz = true);

  site_parameter_hxz() : s_(0.5), c_(0), hx_(0), hz_(0) {}
  site_parameter_hxz(double s) : s_(s), c_(0), hx_(0), hz_(0) {}
  template<typename J>
  site_parameter_hxz(const alps::half_integer<J>& s) :
    s_(s), c_(0), hx_(0), hz_(0) {}
  site_parameter_hxz(double s, double c, double hx, double hz) :
    s_(s), c_(c), hx_(hx), hz_(hz) {}
  template<typename J>
  site_parameter_hxz(const alps::half_integer<J>& s,
                     double c, double hx, double hz) :
    s_(s), c_(c), hx_(hx), hz_(hz) {}
  site_parameter_hxz(const boost::multi_array<double, 2>& mat)
  {
    bool success = fit2site(mat, *this);
    if (!success)
      boost::throw_exception(std::runtime_error("Error: fitting to site_parameter_hxz failed.  This model is not supported by the current looper code."));
  }

  bool operator==(const site_parameter_hxz& rhs) const
  {
    return (s() == rhs.s()) && nearly_equal(c(), rhs.c()) &&
      nearly_equal(hx(), rhs.hx()) && nearly_equal(hz(), rhs.hz());
  }

  bool operator!=(const site_parameter_hxz& rhs) const
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


class site_parameter_noh
{
public:
  typedef alps::half_integer<int> spin_type;

  BOOST_STATIC_CONSTANT(bool, has_hx = false);
  BOOST_STATIC_CONSTANT(bool, has_hz = false);

  site_parameter_noh() : s_(0.5), c_(0) {}
  site_parameter_noh(double s) : s_(s), c_(0) {}
  template<typename J>
  site_parameter_noh(const alps::half_integer<J>& s) :
    s_(s), c_(0) {}
  site_parameter_noh(double s, double c) :
    s_(s), c_(c) {}
  template<typename J>
  site_parameter_noh(const alps::half_integer<J>& s, double c) :
    s_(s), c_(c) {}
  site_parameter_noh(const boost::multi_array<double, 2>& mat)
  {
    site_parameter_hxz p;
    bool success = fit2site(mat, p);
    if (!(success && nearly_equal(p.hx(), 0) && nearly_equal(p.hz(), 0)))
      boost::throw_exception(std::runtime_error("Error: fitting to site_parameter_noh failed.  This model is not supported by the current looper code."));
    s_ = p.s();
    c_ = p.c();
  }

  bool operator==(const site_parameter_noh& rhs) const
  {
    return (s() == rhs.s()) && nearly_equal(c(), rhs.c());
  }

  bool operator!=(const site_parameter_noh& rhs) const
  {
    return !operator==(rhs);
  }

  const spin_type& s() const { return s_; }
  spin_type& s() { return s_; }

  double c() const { return c_; }
  double& c() { return c_; }

private:
  spin_type s_;
  double c_;
};


class bond_parameter_xxz
{
public:
  bond_parameter_xxz() : c_(), jxy_(), jz_() {}
  bond_parameter_xxz(double c, double jxy, double jz) 
    : c_(c), jxy_(jxy), jz_(jz) {}
  bond_parameter_xxz(const boost::multi_array<double, 4>& mat)
  {
    bool success = fit2bond(mat, *this);
    if (!success)
      boost::throw_exception(std::runtime_error("Error: fitting to bond_parameter_xxz failed.  This model is not supported by the corrent looper code."));
  }

  double c() const { return c_; }
  double& c() { return c_; }

  double jxy() const { return jxy_; }
  double& jxy() { return jxy_; }

  double jz() const { return jz_; }
  double& jz() { return jz_; }

  bool operator==(const bond_parameter_xxz& rhs) const
  {
    return nearly_equal(c(), rhs.c()) && nearly_equal(jxy(), rhs.jxy()) &&
      nearly_equal(jz(), rhs.jz());
  }
  bool operator!=(const bond_parameter_xxz& rhs) const
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
  site_matrix(const site_parameter_hxz& s) : matrix_()
  { build(s.s(), s.c(), s.hx(), s.hz()); }

  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  template<typename I>
  void build(const alps::half_integer<I>& s, value_type c, value_type hx,
             value_type hz)
  {
    using std::sqrt; using alps::to_double;

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
      matrix_[sz.distance(-s)][sz.distance(-s)] = c - hz * to_double(sz);

    // off-diagonal elements: hx s+ / 2
    for (half_integer_type sz = -s; sz <= s-1; ++sz)
      matrix_[sz.distance(-s)+1][sz.distance(-s)] =
        - 0.5 * hx * sqrt(to_double(s-sz) * to_double(s+sz+1));

    // off-diagonal elements: hx s- / 2
    for (half_integer_type sz = -s+1; sz <= s; ++sz)
      matrix_[sz.distance(-s)-1][sz.distance(-s)] = 
        - 0.5 * hx * sqrt(to_double(s+sz) * to_double(s-sz+1));
  }

private:
  matrix_type matrix_;
};


template<typename T>
class bond_matrix
{
public:
  typedef T value_type;
  typedef boost::multi_array<value_type, 4> matrix_type;
  typedef typename matrix_type::size_type size_type;

  bond_matrix() : matrix_() {}
  bond_matrix(const bond_matrix& m) : matrix_(m.matrix_) {}
  template<typename I>
  bond_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
              value_type e0, value_type jxy, value_type jz) : matrix_()
  { build(s0, s1, e0, jxy, jz); }
  template<typename SITE_P>
  bond_matrix(const SITE_P& s0, const SITE_P& s1,
             value_type e0, value_type jxy, value_type jz) : matrix_()
  { build(s0.s(), s1.s(), e0, jxy, jz); }
  template<typename SITE_P, typename BOND_P>
  bond_matrix(const SITE_P& s0, const SITE_P& s1, const BOND_P& p) : matrix_()
  { build(s0.s(), s1.s(), p.c(), p.jxy(), p.jz()); }

  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  template<typename I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
             value_type c, value_type jxy, value_type jz)
  {
    using std::sqrt; using alps::to_double;

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
          c - jz * to_double(sz0) * to_double(sz1);
      }
    }

    // off-diagonal elements: - jxy s0+ s1- / 2
    for (half_integer_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
      for (half_integer_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
        matrix_[sz0.distance(-s0)+1][sz1.distance(-s1)-1]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          - 0.5 * jxy *
          sqrt(to_double(s0-sz0) * to_double(s0+sz0+1)) *
          sqrt(to_double(s1+sz1) * to_double(s1-sz1+1));
      }
    }

    // off-diagonal elements: - jxy s0- s1+ / 2
    for (half_integer_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
      for (half_integer_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
        matrix_[sz0.distance(-s0)-1][sz1.distance(-s1)+1]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          - 0.5 * jxy *
          sqrt(to_double(s0+sz0) * to_double(s0-sz0+1)) *
          sqrt(to_double(s1-sz1) * to_double(s1+sz1+1));
      }
    }
  }

private:
  matrix_type matrix_;
};


//
// models
//

template<typename SITE_P, typename BOND_P>
class model_parameter
{
public:
  typedef int type_type;
  typedef SITE_P site_parameter_type;
  typedef BOND_P bond_parameter_type;
  typedef std::map<type_type, site_parameter_type> site_map_type;
  typedef std::map<type_type, bond_parameter_type> bond_map_type;
  typedef typename site_parameter_type::spin_type spin_type;

  template<typename I, typename G>
  model_parameter(double Jxy, double Jz, const alps::half_integer<I>& spin,
                  const G& graph)
    : sites_(), bonds_(), signed_(), frustrated_()
  { set_parameters(Jxy, Jz, spin, graph); }
  template<typename I, typename G>
  model_parameter(double Jxy, double Jz, const alps::half_integer<I>& spin,
                  const G& graph, bool is_signed)
    : sites_(), bonds_(), signed_(), frustrated_()
  { set_parameters(Jxy, Jz, spin, graph, is_signed); }
  template<typename G, typename I>
  model_parameter(const alps::Parameters params, const G& graph,
                  const alps::HamiltonianDescriptor<I>& hd)
    : sites_(), bonds_(), signed_(), frustrated_()
  { set_parameters(params, graph, hd); }
  template<typename G, typename I>
  model_parameter(const alps::Parameters params, const G& graph,
                  const alps::HamiltonianDescriptor<I>& hd,
                  bool is_signed)
    : sites_(), bonds_(), signed_(), frustrated_()
  { set_parameters(params, graph, hd, is_signed); }
  template<typename G>
  model_parameter(const alps::Parameters params, const G& graph,
                  const alps::ModelLibrary& models)
    : sites_(), bonds_(), signed_(), frustrated_()
  { set_parameters(params, graph, models); }
  template<typename G>
  model_parameter(const alps::Parameters params, const G& graph,
                   const alps::ModelLibrary& models, bool is_signed)
    : sites_(), bonds_(), signed_(), frustrated_()
  { set_parameters(params, graph, models, is_signed); }

  // set_parameters

  template<typename I, typename G>
  void set_parameters(double Jxy, double Jz, const alps::half_integer<I>& spin,
                      const G& graph)
  {
    set_parameters_impl(Jxy, Jz, spin, graph);
    signed_ = check_sign(graph);
    frustrated_ = check_classical_frustration(graph);
  }
  template<typename I, typename G>
  void set_parameters(double Jxy, double Jz, const alps::half_integer<I>& spin,
                      const G& graph, bool is_signed)
  {
    set_parameters_impl(Jxy, Jz, spin, graph);
    signed_ = is_signed;
    frustrated_ = check_classical_frustration(graph);
  }
  template<typename G, typename I>
  void set_parameters(alps::Parameters params, const G& graph,
                      const alps::HamiltonianDescriptor<I>& hd)
  {
    set_parameters_impl(params, graph, hd);
    signed_ = check_sign(graph);
    frustrated_ = check_classical_frustration(graph);
  }
  template<typename G, typename I>
  void set_parameters(alps::Parameters params, const G& graph,
                      const alps::HamiltonianDescriptor<I>& hd,
                      bool is_signed)
  {
    set_parameters_impl(params, graph, hd);
    signed_ = is_signed;
    frustrated_ = check_classical_frustration(graph);
  }
  template<typename G>
  void set_parameters(const alps::Parameters params, const G& graph,
                      const alps::ModelLibrary& models)
  {
    // get Hamilton operator from ModelLibrary
    alps::HamiltonianDescriptor<short>
      hd(models.get_hamiltonian(params["MODEL"]), params);
    alps::Parameters p(params);
    p.copy_undefined(hd.default_parameters());
    set_parameters(p, graph, hd);
  }
  template<typename G>
  void set_parameters(const alps::Parameters params,
                      const G& graph,
                      const alps::ModelLibrary& models,
                      bool is_signed)
  {
    // get Hamilton operator from ModelLibrary
    alps::HamiltonianDescriptor<short>
      hd(models.get_hamiltonian(params["MODEL"]), params);
    alps::Parameters p(params);
    p.copy_undefined(hd.default_parameters());
    set_parameters(p, graph, hd, is_signed);
  }

  int num_site_types() const { return sites_.size(); }
  bool is_uniform_site() const { return num_site_types() == 1; }
  site_parameter_type site(type_type t) const
  { return sites_.find(t)->second; }
  site_parameter_type uniform_site() const
  {
#ifndef NDEBUG
    assert(is_uniform_site());
#endif
    return sites_.begin()->second;
  }

  int num_bond_types() const { return bonds_.size(); }
  bool is_uniform_bond() const { return num_bond_types() == 1; }
  bond_parameter_type bond(type_type t) const
    { return bonds_.find(t)->second; }
  bond_parameter_type uniform_bond() const
  {
#ifndef NDEBUG
    assert(is_uniform_bond());
#endif
    return bonds_.begin()->second;
  }

  bool is_signed() const { return signed_; }
  bool is_classically_frustrated() const { return frustrated_; }

protected:
  template<typename I, typename G>
  void set_parameters_impl(double Jxy, double Jz,
                           const alps::half_integer<I>& spin,
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
      if (!sites_.count(t)) sites_[t].s() = spin;
    }

    // set bond parameters
    typename alps::property_map<alps::bond_type_t, graph_type,
                                type_type>::const_type
      bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));

    edge_iterator ei_end = boost::edges(graph).second;
    for (edge_iterator ei = boost::edges(graph).first; ei != ei_end; ++ei) {
      type_type t = bond_type[*ei];
      if (!bonds_.count(t)) bonds_[t] = bond_parameter_xxz(0.0, Jxy, Jz);
    }
  }

  template<typename G, typename I>
  void set_parameters_impl(alps::Parameters params, const G& graph,
    const alps::HamiltonianDescriptor<I>& hd)
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
      if (!sites_.count(t))
        sites_[t] = site_parameter_type(alps::get_matrix(0.,hd.site_term(t),
          hd.basis().site_basis(t), params));
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
        bond_parameter_type p(alps::get_matrix(0., hd.bond_term(bt), 
          hd.basis().site_basis(st0), hd.basis().site_basis(st1), params));
        if (!bonds_.count(bt)) {
          bonds_[bt] = p;
        } else {
          if (bonds_[bt] != p)
            boost::throw_exception(std::runtime_error("inconsistent bond parameter(s)"));
        }
      }
    }
  }

  template<typename G>
  bool check_sign(const G& graph) const
  {
    std::vector<double> w(boost::num_edges(graph));
    typename alps::property_map<alps::bond_type_t, G, type_type>::const_type
      bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));
    typename boost::graph_traits<G>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
      if (bond(bond_type[*ei]).jxy() > 0.0) {
        w[boost::get(boost::edge_index, graph, *ei)] = 1.0;
      } else if (bond(bond_type[*ei]).jxy() < 0.0) {
        w[boost::get(boost::edge_index, graph, *ei)] = -1.0;
      } else {
        w[boost::get(boost::edge_index, graph, *ei)] = 0.0;
      }
    }
    return alps::is_frustrated(graph,
             boost::make_iterator_property_map(w.begin(),
             boost::get(boost::edge_index, graph)));
  }

  template<typename G>
  bool check_classical_frustration(const G& graph) const
  {
    std::vector<double> w(boost::num_edges(graph));
    typename alps::property_map<alps::bond_type_t, G, type_type>::const_type
      bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));
    typename boost::graph_traits<G>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei) {
      if (bond(bond_type[*ei]).jz() > 0.0) {
        w[boost::get(boost::edge_index, graph, *ei)] = 1.0;
      } else if (bond(bond_type[*ei]).jz() < 0.0) {
        w[boost::get(boost::edge_index, graph, *ei)] = -1.0;
      } else {
        w[boost::get(boost::edge_index, graph, *ei)] = 0.0;
      }
    }
    return alps::is_frustrated(graph,
      boost::make_iterator_property_map(w.begin(),
      boost::get(boost::edge_index, graph)));
  }

private:
  site_map_type sites_;
  bond_map_type bonds_;
  bool signed_;
  bool frustrated_;
};


//
// functions for fitting
//

template<typename T>
bool fit2site(const boost::multi_array<T, 2>& mat, site_parameter_hxz& param,
              T tol)
{
  typedef T value_type;

#ifndef NDEBUG
  assert(mat.shape()[0] == mat.shape()[1]);
#endif

  int dim = mat.shape()[0];
  int m = dim * dim;
  int n = 3;
  
  site_parameter_hxz::spin_type s((double)(dim - 1)/2);

  alps::SiteBasisDescriptor<short> basis(spin_basis(s));
  boost::multi_array<value_type, 2> mat_c(
    alps::get_matrix(value_type(),
                     alps::SiteTermDescriptor("1"),
                     basis, alps::Parameters()));
  boost::multi_array<value_type, 2> mat_hx(
    alps::get_matrix(value_type(),
                     alps::SiteTermDescriptor("-(Splus+Sminus)/2"),
                     basis, alps::Parameters()));
  boost::multi_array<value_type, 2> mat_hz(
    alps::get_matrix(value_type(),
                     alps::SiteTermDescriptor("-Sz"),
                     basis, alps::Parameters()));

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

  value_type c = alps::expression::numeric_cast<value_type>(x(0));
  if (std::abs(c) < tol) c = 0;
  value_type hx = alps::expression::numeric_cast<value_type>(x(1));
  if (std::abs(hx) < tol) hx = 0;
  value_type hz = alps::expression::numeric_cast<value_type>(x(2));
  if (std::abs(hz) < tol) hz = 0;

  if (success) param = site_parameter_hxz(s, c, hx, hz);

  return success;
}


template<typename T>
bool fit2site(const site_matrix<T>& mat, site_parameter_hxz& param, T tol)
{
  return fit2site(mat.matrix(), param, tol);
}


template<typename T>
bool fit2bond(const boost::multi_array<T, 4>& mat, bond_parameter_xxz& param,
              T tol)
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

  alps::SiteBasisDescriptor<short> basis0(spin_basis((double)(d0 - 1)/2));
  alps::SiteBasisDescriptor<short> basis1(spin_basis((double)(d1 - 1)/2));

  boost::multi_array<value_type, 4> mat_c(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor("1"),
                     basis0, basis1));
  boost::multi_array<value_type, 4> mat_jxy(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor(
                       "-(Splus(i)*Sminus(j)+Sminus(i)*Splus(j))/2"),
                     basis0, basis1));
  boost::multi_array<value_type, 4> mat_jz(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor("-Sz(i)*Sz(j)"),
                     basis0, basis1));

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

  value_type c = alps::expression::numeric_cast<value_type>(x(0));
  if (std::abs(c) < tol) c = 0;
  value_type jxy = alps::expression::numeric_cast<value_type>(x(1));
  if (std::abs(jxy) < tol) jxy = 0;
  value_type jz = alps::expression::numeric_cast<value_type>(x(2));
  if (std::abs(jz) < tol) jz = 0;

  if (success) param = bond_parameter_xxz(c, jxy, jz);

  return success;
}


template<typename T>
bool fit2bond(const bond_matrix<T>& mat, bond_parameter_xxz& param, T tol)
{
  return fit2bond(mat.matrix(), param, tol);
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
std::ostream& operator<<(std::ostream& os, const looper::bond_matrix<T>& m)
{
  boost::numeric::ublas::matrix<T> mat;
  flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, const looper::site_parameter_hxz& p)
{
  os << "C = " << p.c() << ", Hx = " << p.hx() << ", Hz = " << p.hz();
  return os;
}

inline
std::ostream& operator<<(std::ostream& os, const looper::bond_parameter_xxz& p)
{
  os << "C = " << p.c() << ", Jxy = " << p.jxy() << ", Jz = " << p.jz();
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace looper
#endif

#endif // LOOPER_MODEL_H
