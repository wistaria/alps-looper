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

#include <alps/parameters.h>
#include <alps/lattice.h>
#include <alps/model.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace looper {

template<class I = int>
class site_parameter
{
public:
  typedef alps::half_integer<I> spin_type;

  site_parameter()
    : s_(0.5), c_(), hx_(), hz_() {}
  site_parameter(double s)
    : s_(s), c_(), hx_(), hz_() {}
  template<class J>
  site_parameter(const alps::half_integer<J>& s)
    : s_(s), c_(), hx_(), hz_() {}
  template<class J>
  site_parameter(const alps::half_integer<J>& s,
                 double c, double hx, double hz)
    : s_(s), c_(c), hx_(hx), hz_(hz) {}

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

  double c() const { return c_; }
  double& c() { return c_; }

  double jxy() const { return jxy_; }
  double& jxy() { return jxy_; }

  double jz() const { return jz_; }
  double& jz() { return jz_; }

  bool operator==(const xxz_parameter& rhs) const {
    return (c_ == rhs.c_) && (jxy_ == rhs.jxy_) && (jz_ == rhs.jz_);
  }
  bool operator!=(const xxz_parameter& rhs) const {
    return !operator==(rhs);
  }

private:
  double c_;
  double jxy_;
  double jz_;
};


template<class T, class U>
inline void flatten_matrix(const boost::multi_array<T, 4>& m_in,
                           boost::numeric::ublas::matrix<U>& m_out)
{
#ifndef NDEBUG
  assert(m_in.shape()[0] == m_in.shape()[2]);
  assert(m_in.shape()[1] == m_in.shape()[3]);
#endif

  int d0 = m_in.shape()[0];
  int d1 = m_in.shape()[1];
  int dim = d0 * d1;

  m_out.resize(dim, dim);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
        for (int j1 = 0; j1 < d1; ++j1)
          m_out(i0 * d1 + i1, j0 * d1 + j1) = m_in[i0][i1][j0][j1];
}


template <class T = double>
class xxz_matrix
{
public:
  typedef T value_type;
  typedef boost::multi_array<value_type, 4> matrix_type;
  typedef typename matrix_type::size_type size_type;

  xxz_matrix() : matrix_() {}
  xxz_matrix(const xxz_matrix& m) : matrix_(m.matrix_) {}
  template<class I>
  xxz_matrix(const site_parameter<I>& s0, const site_parameter<I>& s1,
             double e0, double jxy, double jz) : matrix_()
  { build(s0.s(), s1.s(), e0, jxy, jz); }
  template<class I>
  xxz_matrix(const site_parameter<I>& s0, const site_parameter<I>& s1,
             const xxz_parameter& p) : matrix_()
  { build(s0.s(), s1.s(), p.c(), p.jxy(), p.jz()); }

  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  template<class I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
             double e0, double jxy, double jz)
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

    // diagonal elements: jz sz0 sz1
    for (half_integer_type sz0 = s0; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1; --sz1) {
        matrix_[s0.distance(sz0)][s1.distance(sz1)]
          [s0.distance(sz0)][s1.distance(sz1)] =
          e0 - jz * double(sz0) * double(sz1);
      }
    }

    // off-diagonal elements: jxy s0+ s1- / 2
    for (half_integer_type sz0 = s0-1; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1+1; --sz1) {
        matrix_[s0.distance(sz0+1)][s1.distance(sz1-1)]
          [s0.distance(sz0)][s1.distance(sz1)] =
          - 0.5 * jxy *
          std::sqrt(double(s0-sz0) * double(s0+sz0+1)) *
          std::sqrt(double(s1+sz1) * double(s1-sz1+1));
      }
    }

    // off-diagonal elements: jxy s0- s1+ / 2
    for (half_integer_type sz0 = s0; sz0 >= -s0 + 1; --sz0) {
      for (half_integer_type sz1 = s1-1; sz1 >= -s1; --sz1) {
        matrix_[s0.distance(sz0-1)][s1.distance(sz1+1)]
          [s0.distance(sz0)][s1.distance(sz1)] =
          - 0.5 * jxy *
          std::sqrt(double(s0+sz0) * double(s0-sz0+1)) *
          std::sqrt(double(s1-sz1) * double(s1+sz1+1));
      }
    }
  }

private:
  matrix_type matrix_;
};


//
// fitting a matrix to xxz_matrix
//

template <class T>
std::pair<bool, xxz_parameter>
fit2xxz(const boost::multi_array<T, 4>& mat, T tol = 1.0e-10)
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

  alps::SiteBasisDescriptor<short> basis0(spin_basis((double)(d0 - 1)/2.));
  alps::SiteBasisDescriptor<short> basis1(spin_basis((double)(d1 - 1)/2.));

  boost::multi_array<value_type, 4> mat_c(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor<short>("1"),
                     basis0, basis1, spin_operators(), alps::Parameters()));
  boost::multi_array<value_type, 4> mat_jxy(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor<short>(
                       "-(Splus(i)*Sminus(j)+Sminus(i)*Splus(j))/2"),
                     basis0, basis1, spin_operators(), alps::Parameters()));
  boost::multi_array<value_type, 4> mat_jz(
    alps::get_matrix(value_type(),
                     alps::BondTermDescriptor<short>("-Sz(i)*Sz(j)"),
                     basis0, basis1, spin_operators(), alps::Parameters()));

  boost::numeric::ublas::matrix<value_type> a(m, n);
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

  // call linear least sqaure problem solver in LAPACK
  solve_llsp(a, b, x);

  value_type e0 = ((std::abs(x(0)) > tol) ? x(0) : 0.);
  value_type jxy = ((std::abs(x(1)) > tol) ? x(1) : 0.);
  value_type jz = ((std::abs(x(2)) > tol) ? x(2) : 0.);

  // check
  bool success = true;
  for (int i0 = 0; i0 < d0; ++i0) {
    for (int i1 = 0; i1 < d1; ++i1) {
      for (int j0 = 0; j0 < d0; ++j0) {
        for (int j1 = 0; j1 < d1; ++j1) {
          if (std::abs(mat[i0][i1][j0][j1] -
                       (e0 * mat_c[i0][i1][j0][j1] +
                        jxy * mat_jxy[i0][i1][j0][j1] +
                        jz * mat_jz[i0][i1][j0][j1])) > tol) success = false;
        }
      }
    }
  }

  return std::make_pair(success, xxz_parameter(e0, jxy, jz));
}


template <class T>
std::pair<bool, xxz_parameter>
fit2xxz(const xxz_matrix<T>& mat, T tol = 1.0e-10)
{
  return fit2xxz(mat.matrix(), tol);
}


template<class I = int>
class xxz_model
{
public:
  typedef int type_type;
  typedef std::map<type_type, site_parameter<I> > site_type;
  typedef std::map<type_type, xxz_parameter> bond_type;
  typedef typename site_parameter<I>::spin_type spin_type;

  template<class IntType,class G>
  xxz_model(double Jxy, double Jz, const alps::half_integer<IntType>& spin,
            const G& graph) : site_(), bond_()
  { set_parameters(Jxy, Jz, spin, graph); }
  template<class G, class IntType>
  xxz_model(const alps::Parameters params, const G& graph,
            const alps::ModelLibrary::OperatorDescriptorMap& ops,
            const alps::HamiltonianDescriptor<IntType>& hd)
    : site_(), bond_()
  { set_parameters(params, graph, ops, hd); }
  template<class G>
  xxz_model(const alps::Parameters params, const G& graph,
            const alps::ModelLibrary& models) : site_(), bond_()
  { set_parameters(params, graph, models); }

  template<class IntType, class G>
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
      if (!site_.count(t))
        site_[t].s() = (double(hd.basis().site_basis(t).num_states()) - 1) / 2;
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
        bool success;
        xxz_parameter p;
        boost::tie(success, p) =
          fit2xxz(alps::get_matrix(0.,hd.bond_term(bt),hd.basis().site_basis(st0),
                                          hd.basis().site_basis(st1),
                                          ops, params));
        if (!success)
          boost::throw_exception(std::runtime_error("fitting to XXZ model "
                                                    "failed"));
        if (!bond_.count(bt)) {
          bond_[bt] = p;
        } else {
          if (bond_[bt] != p)
            boost::throw_exception(std::runtime_error("inconsistent bond "
                                                      "parameter(s)"));
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


} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline
std::ostream& operator<<(std::ostream& os, const looper::xxz_parameter& p)
{
  os << "C = " << p.c() << ", Jxy = " << p.jxy() << ", Jz = " << p.jz();
  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const looper::xxz_matrix<T>& m)
{
  boost::numeric::ublas::matrix<T> mat;
  flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace looper
#endif

#endif // LOOPER_XXZ_H
