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
  site_parameter_hxz& param, T tol = 1.e-10);
// template<typename T> bool fit2site(const site_matrix<T>& mat,
//   site_parameter_hxz& param, T tol = 1.e-10);
template<typename T> bool fit2bond(const boost::multi_array<T, 4>& mat,
  bond_parameter_xxz& param, T tol = 1.e-10);
// template<typename T> bool fit2bond(const bond_matrix<T>& mat,
//   bond_parameter_xxz& param, T tol = 1.e-10);

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
      boost::throw_exception(std::runtime_error(
        "Error: fitting to site_parameter_hxz failed.  "
	"This model is not supported by the current looper code."));
  }

  bool operator==(const site_parameter_hxz& rhs) const
  {
    return (s() == rhs.s()) && equal(c(), rhs.c()) &&
      equal(hx(), rhs.hx()) && equal(hz(), rhs.hz());
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
    if (!(success && equal(p.hx(), 0.) && equal(p.hz(), 0.)))
      boost::throw_exception(std::runtime_error(
        "Error: fitting to site_parameter_noh failed.  "
	"This model is not supported by the current looper code."));
    s_ = p.s();
    c_ = p.c();
  }

  bool operator==(const site_parameter_noh& rhs) const
  {
    return (s() == rhs.s()) && equal(c(), rhs.c());
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
      boost::throw_exception(std::runtime_error(
        "Error: fitting to bond_parameter_xxz failed.  "
	"This model is not supported by the corrent looper code."));
  }

  double c() const { return c_; }
  double& c() { return c_; }

  double jxy() const { return jxy_; }
  double& jxy() { return jxy_; }

  double jz() const { return jz_; }
  double& jz() { return jz_; }

  bool operator==(const bond_parameter_xxz& rhs) const
  {
    return equal(c(), rhs.c()) && equal(jxy(), rhs.jxy()) &&
      equal(jz(), rhs.jz());
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
  typedef boost::multi_array<T, 2>        matrix_type;
  typedef T                               value_type;
  typedef typename matrix_type::size_type size_type;

  site_matrix() : mat_() {}
  site_matrix(const site_matrix& m) : mat_(m.mat_) {}
  template<typename SITE_P>
  site_matrix(const SITE_P& sp) : mat_()
  {
    using std::sqrt; using alps::to_double;
    typedef typename SITE_P::spin_type spin_type;

    spin_type s = sp.s();

    // set matrix dimension
    int dim = s.get_twice()+1;
    mat_.resize(boost::extents[dim][dim]);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        mat_[i][j] = value_type(0);

    // diagonal elements: c - hz sz
    for (spin_type sz = -s; sz <= s; ++sz)
      mat_[sz.distance(-s)][sz.distance(-s)] =
	sp.c() - sp.hz() * to_double(sz);
    
    // off-diagonal elements: hx s+ / 2
    for (spin_type sz = -s; sz <= s-1; ++sz)
      mat_[sz.distance(-s)+1][sz.distance(-s)] =
	- 0.5 * sp.hx() * sqrt(to_double(s-sz) * to_double(s+sz+1));
    
    // off-diagonal elements: hx s- / 2
    for (spin_type sz = -s+1; sz <= s; ++sz)
      mat_[sz.distance(-s)-1][sz.distance(-s)] = 
	- 0.5 * sp.hx() * sqrt(to_double(s+sz) * to_double(s-sz+1));
  }

  value_type& operator()(size_type i, size_type j) { return mat_[i][j]; }
  const value_type& operator()(size_type i, size_type j) const
  { return mat_[i][j]; }

  matrix_type& matrix() { return mat_; }
  const matrix_type& matrix() const { return mat_; }

private:
  boost::multi_array<T, 2> mat_;
};


template<typename T>
class bond_matrix
{
public:
  typedef boost::multi_array<T, 4>        matrix_type;
  typedef T                               value_type;
  typedef typename matrix_type::size_type size_type;

  bond_matrix() : mat_() {}
  bond_matrix(const bond_matrix& m) : mat_(m.mat_) {}
  template<typename I, typename BOND_P>
  bond_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
              const BOND_P& bp) : mat_()
  { build(s0, s1, bp); }
  template<typename SITE_P, typename BOND_P>
  bond_matrix(const SITE_P& sp0, const SITE_P& sp1, const BOND_P& bp) :
    mat_()
  { build(sp0.s(), sp1.s(), bp); }

  value_type&
  operator()(size_type i, size_type j, size_type k, size_type l)
  { return mat_[i][j][k][l]; }
  const value_type&
  operator()(size_type i, size_type j, size_type k, size_type l) const
  { return mat_[i][j][k][l]; }

  matrix_type& matrix() { return mat_; }
  const matrix_type& matrix() const { return mat_; }

protected:
  template<typename I, typename BOND_P>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
             const BOND_P& bp)
  {
    using std::sqrt; using alps::to_double;
    typedef typename alps::half_integer<I> spin_type;

    // set matrix dimension
    int d0 = s0.get_twice()+1;
    int d1 = s1.get_twice()+1;
    mat_.resize(boost::extents[d0][d1][d0][d1]);
    for (int i0 = 0; i0 < d0; ++i0)
      for (int i1 = 0; i1 < d1; ++i1)
        for (int j0 = 0; j0 < d0; ++j0)
          for (int j1 = 0; j1 < d1; ++j1)
            mat_[i0][i1][j0][j1] = value_type(0);

    // diagonal elements: c - jz sz0 sz1
    for (spin_type sz0 = -s0; sz0 <= s0; ++sz0) {
      for (spin_type sz1 = -s1; sz1 <= s1; ++sz1) {
        mat_[sz0.distance(-s0)][sz1.distance(-s1)]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          bp.c() - bp.jz() * to_double(sz0) * to_double(sz1);
      }
    }

    // off-diagonal elements: - jxy s0+ s1- / 2
    for (spin_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
      for (spin_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
        mat_[sz0.distance(-s0)+1][sz1.distance(-s1)-1]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          - 0.5 * bp.jxy() *
          sqrt(to_double(s0-sz0) * to_double(s0+sz0+1)) *
          sqrt(to_double(s1+sz1) * to_double(s1-sz1+1));
      }
    }

    // off-diagonal elements: - jxy s0- s1+ / 2
    for (spin_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
      for (spin_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
        mat_[sz0.distance(-s0)-1][sz1.distance(-s1)+1]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          - 0.5 * bp.jxy() *
          sqrt(to_double(s0+sz0) * to_double(s0-sz0+1)) *
          sqrt(to_double(s1-sz1) * to_double(s1+sz1+1));
      }
    }
  }

private:
  boost::multi_array<T, 4> mat_;
};


//
// models
//

template<typename SITE_P, typename BOND_P>
class model_parameter
{
public:
  typedef SITE_P site_parameter_type;
  typedef BOND_P bond_parameter_type;
  typedef std::vector<site_parameter_type> site_map_type;
  typedef std::vector<bond_parameter_type> bond_map_type;
  typedef typename site_parameter_type::spin_type spin_type;

  template<typename G, typename I>
  model_parameter(const G& g, const alps::half_integer<I>& spin,
		  double Jxy, double Jz)
    : sites_(), bonds_()
  { set_parameters(g, spin, Jxy, Jz); }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params,
                  const G& g, bool disordered_sites, bool disordered_bond,
		  const alps::model_helper<I>& mh)
    : sites_(), bonds_()
  { set_parameters(params, g, disordered_sites, disordered_bond, mh); }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params,
                  const G& g, bool disordered_sites, bool disordered_bond,
		  const alps::model_helper<I>& mh, bool is_signed)
    : sites_(), bonds_()
  { set_parameters(params, g, disordered_sites, disordered_bond, mh, is_signed); }

  // set_parameters

  template<typename G, typename I>
  void set_parameters(const G& g, const alps::half_integer<I>& spin, 
		      double Jxy, double Jz)
  {
    set_parameters_impl(g, spin, Jxy, Jz);
    signed_ = check_sign(g);
    frustrated_ = check_classical_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params,
                      const G& g, bool disordered_sites, bool disordered_bond,
		      const alps::model_helper<I>& mh)
  {
    set_parameters_impl(params, g, disordered_sites, disordered_bond, mh);
    signed_ = check_sign(g);
    frustrated_ = check_classical_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params,
                      const G& g, bool disordered_sites, bool disordered_bond,
		      const alps::model_helper<I>& mh, bool is_signed)
  {
    set_parameters_impl(params, g, disordered_sites, disordered_bond, mh);
    signed_ = is_signed;
    frustrated_ = check_classical_frustration(g);
  }

  bool uniform_site() const { return sites_.size() == 1; }
  bool disordered_site() const { return use_site_index_; }
  template<class G>
  site_parameter_type site(
    const typename boost::graph_traits<G>::vertex_descriptor& v,
    const G& g) const
  {
    return disordered_site() ?
      sites_[boost::get(vertex_index_t(), g, v)] : 
      (uniform_site() ? sites_[0] :
       sites_[boost::get(vertex_type_t(), g, v)]);
  }

  bool uniform_bond() const { return bonds_.size() == 1; }
  bool disordered_bond() const { return use_bond_index_; }
  template<class G>
  bond_parameter_type bond(
    const typename alps::graph_traits<G>::edge_descriptor& e,
    const G& g) const
  {
    return disordered_bond() ?
      bonds_[boost::get(edge_index_t(), g, e)] :
      (uniform_bond() ? bonds_[0] :
       bonds_[boost::get(edge_type_t(), g, e)]);
  }

  bool is_signed() const { return signed_; }
  bool is_classically_frustrated() const { return frustrated_; }

protected:
  template<typename G, typename I>
  void set_parameters_impl(const G& /* g */, const alps::half_integer<I>& spin,
			   double Jxy, double Jz)
  {
    // set site parameters
    use_site_index_ = false;
    sites_.resize(1);
    sites_[0].s() = spin;

    // set bond parameters
    use_bond_index_ = false;
    bonds_.resize(1);
    bonds_[0] = bond_parameter_xxz(0., Jxy, Jz);
  }

  template<typename G, typename I>
  void set_parameters_impl(alps::Parameters params, const G& g,
    bool disordered_sites, bool disordered_bond,
    const alps::model_helper<I>& mh)
  {
    using boost::get;
    
    typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<G>::edge_iterator edge_iterator;

    params.copy_undefined(mh.model().default_parameters());
    alps::basis_states_descriptor<I> basis(mh.model().basis(), g);
    alps::Disorder::seed(params.value_or_default("DISORDER_SEED",0));

    //
    // site terms
    //
    
    // check type range and resize 'sites_'
    bool use_site_index_ = false;
    if (disordered_sites) {
      use_site_index_ = true;
      sites_.resize(boost::num_vertices(g));
    } else {
      alps::type_type type_min = 0;
      alps::type_type type_max = 0;
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
        type_min = std::min(type_min, get(vertex_type_t(), g, *vi));
        type_max = std::max(type_min, get(vertex_type_t(), g, *vi));
      }
      if (is_negative(type_min) || type_max >= boost::num_vertices(g)) {
        use_site_index_ = true;
        sites_.resize(boost::num_vertices(g));
      } else {
        sites_.resize(type_max);
      }
    }

    // generate site matrices and set site parameters
    if (use_site_index_) {
      alps::Parameters p(params);
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
        if (disordered_sites) {
          alps::throw_if_xyz_defined(params, g);
          p << alps::coordinate_as_parameter(g, *vi);
        }
	unsigned int i = get(vertex_index_t(), g, *vi);
	unsigned int t = get(vertex_type_t(), g, *vi);
        sites_[i] = site_parameter_type(
          alps::get_matrix(double(), mh.model().site_term(t),
                           mh.model().basis().site_basis(t), p));
      }
    } else {
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
	unsigned int t = get(vertex_type_t(), g, *vi);
        sites_[t] = site_parameter_type(
          alps::get_matrix(double(), mh.model().site_term(t),
                           mh.model().basis().site_basis(t), params));
      }
    }

    //
    // bond terms
    //

    // check type range and resize 'bonds_'
    use_bond_index_ = false;
    if (disordered_bond) {
      use_bond_index_ = true;
      bonds_.resize(boost::num_edges(g));
    } else {
      alps::type_type type_min = 0;
      alps::type_type type_max = 0;
      std::map<alps::type_type, std::pair<alps::type_type, alps::type_type> >
	vtype;
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        type_min = std::min(type_min, get(edge_type_t(), g, *ei));
        type_max = std::max(type_max, get(edge_type_t(), g, *ei));
	unsigned int t = get(edge_type_t(), g, *ei);
        if (vtype.find(t) == vtype.end()) {
          vtype[t] = std::make_pair(
              get(vertex_type_t(), g, boost::source(*ei, g)),
              get(vertex_type_t(), g, boost::target(*ei, g)));
        } else {
          if (vtype[t] != std::make_pair(
              get(vertex_type_t(), g, boost::source(*ei, g)),
              get(vertex_type_t(), g, boost::target(*ei, g)))) {
            use_bond_index_ = true;
            break;
          }
        }
      }
      if (use_bond_index_ || is_negative(type_min) ||
          type_max >= boost::num_edges(g)) {
        use_bond_index_ = true;
        bonds_.resize(boost::num_edges(g));
      } else {
        bonds_.resize(type_max);
      }
    }

    // generate bond matrices and set bond parameters
    if (use_bond_index_) {
      alps::Parameters p(params);
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        if (disordered_bond) {
          alps::throw_if_xyz_defined(params, g);
          p << alps::coordinate_as_parameter(g, *ei);
        }
	unsigned int i = get(edge_index_t(), g, *ei);
	unsigned int t = get(edge_type_t(), g, *ei);
	unsigned int st0 = get(vertex_type_t(), g, boost::source(*ei, g));
	unsigned int st1 = get(vertex_type_t(), g, boost::target(*ei, g));
        bonds_[i] = bond_parameter_type(
          alps::get_matrix(double(), mh.model().bond_term(t),
			   mh.model().basis().site_basis(st0),
			   mh.model().basis().site_basis(st1), p));
      }
    } else {
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
	unsigned int t = get(edge_type_t(), g, *ei);
	unsigned int st0 = get(vertex_type_t(), g, boost::source(*ei, g));
	unsigned int st1 = get(vertex_type_t(), g, boost::target(*ei, g));
        bonds_[t] = bond_parameter_type(
          alps::get_matrix(double(), mh.model().bond_term(t),
			   mh.model().basis().site_basis(st0),
			   mh.model().basis().site_basis(st1), params));
      }
    }
  }

  template<typename G>
  bool check_sign(const G& g) const
  {
    std::vector<double> w(boost::num_edges(g));
    typename boost::graph_traits<G>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
      unsigned int i = boost::get(edge_index_t(), g, *ei);
      if (bond(*ei, g).jxy() > 0.) {
        w[i] = 1.;
      } else if (bond(*ei, g).jxy() < 0.) {
        w[i] = -1.;
      } else {
        w[i] = 0.;
      }
    }
    return alps::is_frustrated(g,
             boost::make_iterator_property_map(w.begin(),
	       boost::get(edge_index_t(), g)));
  }

  template<typename G>
  bool check_classical_frustration(const G& g) const
  {
    std::vector<double> w(boost::num_edges(g));
    typename boost::graph_traits<G>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
      unsigned int i = boost::get(edge_index_t(), g, *ei);
      if (bond(*ei, g).jz() > 0.) {
        w[i] = 1.;
      } else if (bond(*ei, g).jz() < 0.) {
        w[i] = -1.;
      } else {
        w[i] = 0.;
      }
    }
    return alps::is_frustrated(g,
      boost::make_iterator_property_map(w.begin(),
        boost::get(edge_index_t(), g)));
  }

private:
  site_map_type sites_;
  bond_map_type bonds_;
  bool use_site_index_;
  bool use_bond_index_;
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
  assert(mat.shape()[0] == mat.shape()[1]);

  typedef T value_type;

  int dim = mat.shape()[0];
  int m = dim * dim;
  int n = 3;
  
  alps::half_integer<short> s((double)(dim-1)/2);
  site_matrix<> mat_c(site_parameter_hxz(s, 1, 0, 0));
  site_matrix<> mat_hx(site_parameter_hxz(s, 0, 1, 0));
  site_matrix<> mat_hz(site_parameter_hxz(s, 0, 0, 1));

//   alps::SiteBasisDescriptor<short> basis(spin_basis(s));
//   boost::multi_array<value_type, 2> mat_c(
//     alps::get_matrix(value_type(),
//                      alps::SiteTermDescriptor("1"),
//                      basis, alps::Parameters()));
//   boost::multi_array<value_type, 2> mat_hx(
//     alps::get_matrix(value_type(),
//                      alps::SiteTermDescriptor("-(Splus(i)+Sminus(i))/2", "i"),
//                      basis, alps::Parameters()));
//   boost::multi_array<value_type, 2> mat_hz(
//     alps::get_matrix(value_type(),
//                      alps::SiteTermDescriptor("-Sz(i)", "i"),
//                      basis, alps::Parameters()));

  boost::numeric::ublas::matrix<value_type,
    boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<value_type> b(m);
  boost::numeric::ublas::vector<value_type> x(n);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      int k = dim * i + j;
      b(k) = mat[i][j];
      a(k, 0) = mat_c(i,j);
      a(k, 1) = mat_hx(i,j);
      a(k, 2) = mat_hz(i,j);
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

// template<typename T>
// bool fit2site(const site_matrix<T>& mat, site_parameter_hxz& param, T tol)
// {
//   return fit2site(mat.matrix(), param, tol);
// }


template<typename T>
bool fit2bond(const boost::multi_array<T, 4>& mat, bond_parameter_xxz& param,
              T tol)
{
  assert(mat.shape()[0] == mat.shape()[2]);
  assert(mat.shape()[1] == mat.shape()[3]);

  typedef T value_type;

  int d0 = mat.shape()[0];
  int d1 = mat.shape()[1];
  int dim = d0 * d1;
  int m = dim * dim;
  int n = 3;

  alps::half_integer<short> s0((double)(d0-1)/2);
  alps::half_integer<short> s1((double)(d1-1)/2);
  bond_matrix<> mat_c(s0, s1, bond_parameter_xxz(1, 0, 0));
  bond_matrix<> mat_jxy(s0, s1, bond_parameter_xxz(0, 1, 0));
  bond_matrix<> mat_jz(s0, s1, bond_parameter_xxz(0, 0, 1));

//   alps::SiteBasisDescriptor<short> basis0(spin_basis((double)(d0 - 1)/2));
//   alps::SiteBasisDescriptor<short> basis1(spin_basis((double)(d1 - 1)/2));
//   boost::multi_array<value_type, 4> mat_c(
//     alps::get_matrix(value_type(),
//                      alps::BondTermDescriptor("1"),
//                      basis0, basis1));
//   boost::multi_array<value_type, 4> mat_jxy(
//     alps::get_matrix(value_type(),
//                      alps::BondTermDescriptor(
//                        "-(Splus(i)*Sminus(j)+Sminus(i)*Splus(j))/2"),
//                      basis0, basis1));
//   boost::multi_array<value_type, 4> mat_jz(
//     alps::get_matrix(value_type(),
//                      alps::BondTermDescriptor("-Sz(i)*Sz(j)"),
//                      basis0, basis1));

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
          a(k, 0) = mat_c(i0,i1,j0,j1);
          a(k, 1) = mat_jxy(i0,i1,j0,j1);
          a(k, 2) = mat_jz(i0,i1,j0,j1);
        }
      }
    }
  }

  // call linear least squaure problem solver
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


// template<typename T>
// bool fit2bond(const bond_matrix<T>& mat, bond_parameter_xxz& param, T tol)
// {
//   return fit2bond(mat.matrix(), param, tol);
// }


} // end namespace looper


#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<typename T>
std::ostream& operator<<(std::ostream& os, const looper::site_matrix<T>& m)
{
  boost::numeric::ublas::matrix<T> mat;
  flatten_matrix(m.matrix(), mat);
  os << mat;
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
