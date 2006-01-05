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

#include "lapack.h"
#include "lattice.h"
#include "util.h"

#include <alps/math.hpp>
#include <alps/model.h>
#include <alps/scheduler/montecarlo.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace looper {

//
// parameters
//

struct site_parameter
{
  typedef alps::half_integer<int> spin_type;

  spin_type s; // size of spin
  double c;    // constant energy offset: + c
  double hx;   // transverse field:       - hx * Sx
  double hz;   // longitudinal field:     - hz * Sz
  double d;    // single-ion anisotropy:  + d * (Sz)^2

  site_parameter(double s_in = 0.5, double c_in = 0, double hx_in = 0,
                 double hz_in = 0, double d_in = 0) :
    s(s_in), c(c_in), hx(hx_in), hz(hz_in), d(d_in) {}
  template<typename J>
  site_parameter(const alps::half_integer<J>& s_in, double c_in = 0,
                 double hx_in = 0, double hz_in = 0, double d_in = 0) :
    s(s_in), c(c_in), hx(hx_in), hz(hz_in), d(d_in) {}
  site_parameter(const boost::multi_array<double, 2>& mat) { do_fit(mat); }

  bool operator==(const site_parameter& rhs) const
  {
    return (s == rhs.s) && alps::is_equal<1>(c, rhs.c) &&
      alps::is_equal<1>(hx, rhs.hx) && alps::is_equal<1>(hz, rhs.hz) &&
      alps::is_equal<1>(d, rhs.d);
  }

  bool operator!=(const site_parameter& rhs) const
  { return !operator==(rhs); }

private:
  void do_fit(const boost::multi_array<double, 2>& mat);
};


struct bond_parameter
{
  double c;   // constant energy offset:   + c
  double jxy; // off-diagonal interaction: + jxy/2 * (s1+ s2- + hc)
  double jz;  // diagonal interaction:     + jz s1z s2z

  bond_parameter(double c_in = 0, double jxy_in = 0, double jz_in = 0)
    : c(c_in), jxy(jxy_in), jz(jz_in) {}
  bond_parameter(const boost::multi_array<double, 4>& mat) { do_fit(mat); }

  bool operator==(const bond_parameter& rhs) const
  {
    return alps::is_equal<1>(c, rhs.c) &&
      alps::is_equal<1>(jxy, rhs.jxy) &&
      alps::is_equal<1>(jz, rhs.jz);
  }
  bool operator!=(const bond_parameter& rhs) const
  { return !operator==(rhs); }

private:
  void do_fit(const boost::multi_array<double, 4>& mat);
};


//
// matrices
//

class site_matrix
{
public:
  typedef boost::multi_array<double, 2>   matrix_type;
  typedef matrix_type::size_type size_type;

  site_matrix() : mat_() {}
  site_matrix(const site_matrix& m) : mat_(m.mat_) {}
  site_matrix(const site_parameter& sp) : mat_()
  {
    using std::sqrt; // using alps::to_double;
    typedef site_parameter::spin_type spin_type;

    spin_type s = sp.s;

    // set matrix dimension
    int dim = s.get_twice()+1;
    mat_.resize(boost::extents[dim][dim]);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        mat_[i][j] = 0;

    // diagonal elements: c - hz sz + d sz^2
    for (spin_type sz = -s; sz <= s; ++sz)
      mat_[sz.distance(-s)][sz.distance(-s)] =
        sp.c - sp.hz * to_double(sz) +  sp.d * power2(to_double(sz));

    // off-diagonal elements: - hx s+ / 2
    for (spin_type sz = -s; sz <= s-1; ++sz)
      mat_[sz.distance(-s)+1][sz.distance(-s)] =
        - 0.5 * sp.hx * sqrt(to_double(s-sz) * to_double(s+sz+1));

    // off-diagonal elements: - hx s- / 2
    for (spin_type sz = -s+1; sz <= s; ++sz)
      mat_[sz.distance(-s)-1][sz.distance(-s)] =
        - 0.5 * sp.hx * sqrt(to_double(s+sz) * to_double(s-sz+1));
  }

  double operator()(size_type i, size_type j) const { return mat_[i][j]; }

  const matrix_type& matrix() const { return mat_; }

private:
  boost::multi_array<double, 2> mat_;
};


class bond_matrix
{
public:
  typedef boost::multi_array<double, 4>   matrix_type;
  typedef matrix_type::size_type size_type;

  bond_matrix() : mat_() {}
  bond_matrix(const bond_matrix& m) : mat_(m.mat_) {}
  template<typename I>
  bond_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
              const bond_parameter& bp) : mat_()
  { build(s0, s1, bp); }
  bond_matrix(const site_parameter& sp0, const site_parameter& sp1,
              const bond_parameter& bp) : mat_()
  { build(sp0.s, sp1.s, bp); }

  double operator()(size_type i, size_type j, size_type k, size_type l) const
  { return mat_[i][j][k][l]; }

  const matrix_type& matrix() const { return mat_; }

protected:
  template<typename I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
             const bond_parameter& bp)
  {
    using std::sqrt; // using alps::to_double;
    typedef typename alps::half_integer<I> spin_type;

    // set matrix dimension
    int d0 = s0.get_twice()+1;
    int d1 = s1.get_twice()+1;
    mat_.resize(boost::extents[d0][d1][d0][d1]);
    for (int i0 = 0; i0 < d0; ++i0)
      for (int i1 = 0; i1 < d1; ++i1)
        for (int j0 = 0; j0 < d0; ++j0)
          for (int j1 = 0; j1 < d1; ++j1)
            mat_[i0][i1][j0][j1] = 0;

    // diagonal elements: c + jz sz0 sz1
    for (spin_type sz0 = -s0; sz0 <= s0; ++sz0) {
      for (spin_type sz1 = -s1; sz1 <= s1; ++sz1) {
        mat_[sz0.distance(-s0)][sz1.distance(-s1)]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          bp.c + bp.jz * to_double(sz0) * to_double(sz1);
      }
    }

    // off-diagonal elements: jxy s0+ s1- / 2
    for (spin_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
      for (spin_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
        mat_[sz0.distance(-s0)+1][sz1.distance(-s1)-1]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          0.5 * bp.jxy *
          sqrt(to_double(s0-sz0) * to_double(s0+sz0+1)) *
          sqrt(to_double(s1+sz1) * to_double(s1-sz1+1));
      }
    }

    // off-diagonal elements: jxy s0- s1+ / 2
    for (spin_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
      for (spin_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
        mat_[sz0.distance(-s0)-1][sz1.distance(-s1)+1]
          [sz0.distance(-s0)][sz1.distance(-s1)] =
          0.5 * bp.jxy *
          sqrt(to_double(s0+sz0) * to_double(s0-sz0+1)) *
          sqrt(to_double(s1-sz1) * to_double(s1+sz1+1));
      }
    }
  }

private:
  boost::multi_array<double, 4> mat_;
};


//
// models
//

class model_parameter
{
public:
  typedef std::vector<site_parameter> site_map_type;
  typedef std::vector<bond_parameter> bond_map_type;
  typedef site_parameter::spin_type spin_type;

  model_parameter() {}
  template<typename G, typename I>
  model_parameter(const G& g, const alps::half_integer<I>& spin,
                  double Jxy, double Jz) { set_parameters(g, spin, Jxy, Jz); }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params, const G& g,
                  bool inhomogeneous_sites, bool inhomogeneous_bond,
                  const alps::HamiltonianDescriptor<I>& hd)
  { set_parameters(params, g, inhomogeneous_sites, inhomogeneous_bond, hd); }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params, const G& g,
                  bool inhomogeneous_sites, bool inhomogeneous_bond,
                  const alps::HamiltonianDescriptor<I>& hd,
                  bool is_signed)
  {
    set_parameters(params, g, inhomogeneous_sites, inhomogeneous_bond,
                   hd, is_signed);
  }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params,
                  const alps::graph_helper<G>& gh,
                  const alps::model_helper<I>& mh)
  {
    set_parameters(params, gh.graph(), gh.inhomogeneous_sites(),
                   gh.inhomogeneous_bonds(), mh.model(),
                   alps::has_sign_problem(mh.model(), gh, params));
  }
  template<typename G, typename I>
  model_parameter(const alps::Parameters& params,
                  const alps::scheduler::LatticeModelMCRun<G, I>& mcrun)
  { set_parameters(params, mcrun); }

  // set_parameters

  template<typename G, typename I>
  void set_parameters(const G& g, const alps::half_integer<I>& spin,
                      double Jxy, double Jz)
  {
    set_parameters_impl(g, spin, Jxy, Jz);
    signed_ = check_sign(g);
    frustrated_ = check_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params,
                      const G& g,
                      bool inhomogeneous_sites,
                      bool inhomogeneous_bond,
                      const alps::HamiltonianDescriptor<I>& hd)
  {
    set_parameters_impl(params, g, inhomogeneous_sites, inhomogeneous_bond,
                        hd);
    signed_ = check_sign(g);
    frustrated_ = check_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params,
                      const G& g,
                      bool inhomogeneous_sites,
                      bool inhomogeneous_bond,
                      const alps::HamiltonianDescriptor<I>& hd,
                      bool is_signed)
  {
    set_parameters_impl(params, g, inhomogeneous_sites, inhomogeneous_bond,
                        hd);
    signed_ = is_signed;
    frustrated_ = check_frustration(g);
  }
  template<typename G, typename I>
  void set_parameters(const alps::Parameters& params,
                      const alps::scheduler::LatticeModelMCRun<G, I>& mcrun)
  {
    set_parameters(params, mcrun.graph(), mcrun.inhomogeneous_sites(),
                   mcrun.inhomogeneous_bonds(), mcrun.model(),
                   alps::has_sign_problem(mcrun.model(), mcrun, params));
  }


  bool uniform_site() const { return sites_.size() == 1; }
  bool inhomogeneous_site() const { return use_site_index_; }
  template<class G>
  const site_parameter& site(
    const typename alps::graph_traits<G>::site_descriptor& v,
    const G& g) const
  {
    return inhomogeneous_site() ?
      sites_[get(site_index_t(), g, v)] :
      (uniform_site() ? sites_[0] :
       sites_[get(site_type_t(), g, v)]);
  }
  const site_parameter& site() const
  {
    if (!uniform_site())
      boost::throw_exception(std::runtime_error("nonuniform sites"));
    return sites_[0];
  }

  bool uniform_bond() const { return bonds_.size() == 1; }
  bool inhomogeneous_bond() const { return use_bond_index_; }
  template<class G>
  const bond_parameter& bond(
    const typename alps::graph_traits<G>::bond_descriptor& e,
    const G& g) const
  {
    return inhomogeneous_bond() ?
      bonds_[get(bond_index_t(), g, e)] :
      (uniform_bond() ? bonds_[0] :
       bonds_[get(bond_type_t(), g, e)]);
  }
  const bond_parameter& bond() const
  {
    if (!uniform_bond())
      boost::throw_exception(std::runtime_error("nonuniform bonds"));
    return bonds_[0];
  }

  bool is_signed() const { return signed_; }
  bool is_frustrated() const { return frustrated_; }
  bool has_field() const { return has_hz_; }
  bool has_d_term() const { return has_d_; }
  double energy_offset() const { return energy_offset_; }

protected:
  template<typename G, typename I>
  void set_parameters_impl(const G& /* g */, const alps::half_integer<I>& spin,
    double Jxy, double Jz);

  template<typename G, typename I>
  void set_parameters_impl(alps::Parameters params, const G& g,
    bool inhomogeneous_sites, bool inhomogeneous_bond,
    const alps::HamiltonianDescriptor<I>& hd);

  template<typename G>
  bool check_sign(const G& g) const;

  template<typename G>
  bool check_frustration(const G& g) const;

private:
  site_map_type sites_;
  bond_map_type bonds_;
  bool use_site_index_;
  bool use_bond_index_;
  bool has_hz_;
  bool has_d_;
  bool signed_;
  bool frustrated_;
  double energy_offset_;
};


////////////////////////////////////////////////////////////////////////////////
//
// Implementations
//

//
// class site_paramter
//

inline void site_parameter::do_fit(const boost::multi_array<double, 2>& mat)
{
  int dim = mat.shape()[0];
  int m = dim * dim;
  int n = (dim > 2 ? 4 : 3); // D-term is meaningful only for S > 1/2
  s = (double)(dim-1)/2;

  if (!(mat.shape()[0] == mat.shape()[1] && s.get_twice() + 1 == dim))
    boost::throw_exception(std::runtime_error(
      "Error: fitting to site_parameter failed.  "
      "This model is not supported by the current version of looper."));

  site_matrix mat_c(site_parameter(s, 1, 0, 0, 0));
  site_matrix mat_hx(site_parameter(s, 0, 1, 0, 0));
  site_matrix mat_hz(site_parameter(s, 0, 0, 1, 0));
  site_matrix mat_d(site_parameter(s, 0, 0, 0, 1));

  boost::numeric::ublas::matrix<double,
     boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<double> b(m);
  boost::numeric::ublas::vector<double> x(n);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      int k = dim * i + j;
      b(k) = mat[i][j];
      a(k, 0) = mat_c(i,j);
      a(k, 1) = mat_hx(i,j);
      a(k, 2) = mat_hz(i,j);
      if (n == 4) a(k, 3) = mat_d(i,j);
    }
  }

  // call linear least sqaure problem solver
  bool success = alps::is_zero<1>(solve_llsp(a, b, x));
  if (!success)
    boost::throw_exception(std::runtime_error(
      "Error: fitting to site_parameter failed.  "
      "This model is not supported by the current version of looper."));

  for (int i = 0; i < n; ++i) x(i) = alps::round<1>(x(i));
  c = x(0);
  hx = x(1);
  hz = x(2);
  d = (n == 4 ? x(3) : 0);
}


//
// class bond_parameter
//

inline void bond_parameter::do_fit(const boost::multi_array<double, 4>& mat)
{
  int d0 = mat.shape()[0];
  int d1 = mat.shape()[1];
  int dim = d0 * d1;
  int m = dim * dim;
  int n = 3;
  alps::half_integer<int> s0((double)(d0-1)/2);
  alps::half_integer<int> s1((double)(d1-1)/2);

  if (!(mat.shape()[0] == mat.shape()[2] && mat.shape()[1] == mat.shape()[3] &&
        s0.get_twice() + 1 == d0 && s1.get_twice() + 1 == d1))
    boost::throw_exception(std::runtime_error(
      "Error: fitting to bond_parameter failed.  "
      "This model is not supported by the current version of looper."));

  bond_matrix mat_c(s0, s1, bond_parameter(1, 0, 0));
  bond_matrix mat_jxy(s0, s1, bond_parameter(0, 1, 0));
  bond_matrix mat_jz(s0, s1, bond_parameter(0, 0, 1));

  boost::numeric::ublas::matrix<double,
    boost::numeric::ublas::column_major> a(m, n);
  boost::numeric::ublas::vector<double> b(m);
  boost::numeric::ublas::vector<double> x(n);
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
  bool success = alps::is_zero<1>(solve_llsp(a, b, x));
  if (!success)
    boost::throw_exception(std::runtime_error(
      "Error: fitting to bond_parameter failed.  "
      "This model is not supported by the current version of looper."));

  for (int i = 0; i < n; ++i) x(i) = alps::round<1>(x(i));
  c = x(0);
  jxy = x(1);
  jz = x(2);
}


//
// class model_parameter
//

template<typename G, typename I>
void model_parameter::set_parameters_impl(const G& /* g */,
  const alps::half_integer<I>& spin, double Jxy, double Jz)
{
  // set site parameters
  use_site_index_ = false;
  sites_.resize(1);
  sites_[0].s = spin;

  // set bond parameters
  use_bond_index_ = false;
  bonds_.resize(1);
  bonds_[0] = bond_parameter(0., Jxy, Jz);

  has_hz_ = false;
  has_d_ = false;
}

template<typename G, typename I>
void model_parameter::set_parameters_impl(alps::Parameters params, const G& g,
  bool inhomogeneous_sites, bool inhomogeneous_bond,
  const alps::HamiltonianDescriptor<I>& hd)
{
  typedef typename alps::graph_traits<G>::site_iterator site_iterator;
  typedef typename alps::graph_traits<G>::bond_iterator bond_iterator;

  params.copy_undefined(hd.default_parameters());
  alps::basis_states_descriptor<I> basis(hd.basis(), g);
  alps::Disorder::seed(params.value_or_default("DISORDER_SEED",0));

  // site terms

  // check type range and resize 'sites_'
  use_site_index_ = false;
  if (inhomogeneous_sites) {
    use_site_index_ = true;
    sites_.resize(num_sites(g));
  } else {
    alps::type_type type_min = 0;
    alps::type_type type_max = 0;
    site_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = sites(g); vi != vi_end; ++vi) {
      type_min = std::min(type_min, get(site_type_t(), g, *vi));
      type_max = std::max(type_min, get(site_type_t(), g, *vi));
    }
    if (alps::is_negative(type_min) || type_max > num_sites(g)) {
      use_site_index_ = true;
      sites_.resize(num_sites(g));
    } else {
      sites_.resize(type_max+1);
    }
  }

  // generate site matrices and set site parameters
  if (use_site_index_) {
    alps::Parameters p(params);
    site_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = sites(g); vi != vi_end; ++vi) {
      if (inhomogeneous_sites) {
        alps::throw_if_xyz_defined(params, g);
        p << alps::coordinate_as_parameter(g, *vi);
      }
      unsigned int i = get(site_index_t(), g, *vi);
      unsigned int t = get(site_type_t(), g, *vi);
      sites_[i] = site_parameter(get_matrix(double(), hd.site_term(t),
        hd.basis().site_basis(t), p));
      if (alps::is_nonzero<1>(sites_[i].hz)) has_hz_ = true;
      if (alps::is_nonzero<1>(sites_[i].d)) has_d_ = true;
    }
  } else {
    std::vector<bool> checked(sites_.size(), false);
    site_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = sites(g); vi != vi_end; ++vi) {
      unsigned int t = get(site_type_t(), g, *vi);
      if (!checked[t]) {
        sites_[t] = site_parameter(get_matrix(double(), hd.site_term(t),
          hd.basis().site_basis(t), params));
        if (alps::is_nonzero<1>(sites_[t].hz)) has_hz_ = true;
        if (alps::is_nonzero<1>(sites_[t].d)) has_d_ = true;
        checked[t] = true;
      }
    }
  }

  // bond terms

  // check type range and resize 'bonds_'
  use_bond_index_ = false;
  if (inhomogeneous_bond) {
    use_bond_index_ = true;
    bonds_.resize(num_bonds(g));
  } else {
    alps::type_type type_min = 0;
    alps::type_type type_max = 0;
    std::map<alps::type_type, std::pair<alps::type_type, alps::type_type> >
      vtype;
    bond_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = bonds(g); ei != ei_end; ++ei) {
      type_min = std::min(type_min, get(bond_type_t(), g, *ei));
      type_max = std::max(type_max, get(bond_type_t(), g, *ei));
      unsigned int t = get(bond_type_t(), g, *ei);
      if (vtype.find(t) == vtype.end()) {
        vtype[t] = std::make_pair(
          get(site_type_t(), g, source(*ei, g)),
          get(site_type_t(), g, target(*ei, g)));
      } else {
        if (vtype[t] != std::make_pair(
          get(site_type_t(), g, source(*ei, g)),
          get(site_type_t(), g, target(*ei, g)))) {
          use_bond_index_ = true;
          break;
        }
      }
    }
    if (use_bond_index_ || alps::is_negative(type_min) ||
        type_max > num_bonds(g)) {
      use_bond_index_ = true;
      bonds_.resize(num_bonds(g));
    } else {
      bonds_.resize(type_max+1);
    }
  }

  // generate bond matrices and set bond parameters
  if (use_bond_index_) {
    alps::Parameters p(params);
    bond_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = bonds(g); ei != ei_end; ++ei) {
      if (inhomogeneous_bond) {
        alps::throw_if_xyz_defined(params, g);
        p << alps::coordinate_as_parameter(g, *ei);
      }
      unsigned int i = get(bond_index_t(), g, *ei);
      unsigned int t = get(bond_type_t(), g, *ei);
      unsigned int st0 = get(site_type_t(), g, source(*ei, g));
      unsigned int st1 = get(site_type_t(), g, target(*ei, g));
      bonds_[i] = bond_parameter(get_matrix(double(), hd.bond_term(t),
        hd.basis().site_basis(st0),
        hd.basis().site_basis(st1), p));
    }
  } else {
    std::vector<bool> checked(bonds_.size(), false);
    bond_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = bonds(g); ei != ei_end; ++ei) {
      unsigned int t = get(bond_type_t(), g, *ei);
      if (!checked[t]) {
        unsigned int st0 = get(site_type_t(), g, source(*ei, g));
        unsigned int st1 = get(site_type_t(), g, target(*ei, g));
        bonds_[t] = bond_parameter(get_matrix(double(), hd.bond_term(t),
          hd.basis().site_basis(st0),
          hd.basis().site_basis(st1), params));
        checked[t] = true;
      }
    }
  }

  energy_offset_ = 0;
  for (site_iterator si = sites(g).first; si != sites(g).second; ++si)
    energy_offset_ += site(*si, g).c;
  for (bond_iterator bi = bonds(g).first; bi != bonds(g).second; ++bi)
    energy_offset_ += bond(*bi, g).c;
  if (has_d_term())
    for (site_iterator si = sites(g).first; si != sites(g).second; ++si)
      energy_offset_ += 0.5 * site(*si, g).s.get_twice() * site(*si, g).d;
}

template<typename G>
bool model_parameter::check_sign(const G& g) const
{
  std::vector<double> w(num_bonds(g));
  typename alps::graph_traits<G>::bond_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = bonds(g); ei != ei_end; ++ei) {
    unsigned int i = get(bond_index_t(), g, *ei);
    if (alps::is_zero<1>(bond(*ei, g).jxy)) {
      w[i] = 0;
    } else if (bond(*ei, g).jxy > 0) {
      w[i] = 1;
    } else {
      w[i] = -1;
    }
  }
  return alps::is_frustrated(g, boost::make_iterator_property_map(w.begin(),
    get(bond_index_t(), g)));
}

template<typename G>
bool model_parameter::check_frustration(const G& g) const
{
  std::vector<double> w(num_bonds(g));
  typename alps::graph_traits<G>::bond_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = bonds(g); ei != ei_end; ++ei) {
    unsigned int i = get(bond_index_t(), g, *ei);
    if (alps::is_zero<1>(bond(*ei, g).jz)) {
      w[i] = 0;
    } else if (bond(*ei, g).jz > 0) {
      w[i] = 1;
    } else {
      w[i] = -1;
    }
  }
  bool frustrated =
    alps::is_frustrated(g, boost::make_iterator_property_map(w.begin(),
      get(bond_index_t(), g)));
  if (has_d_ && !frustrated) {
    typename alps::graph_traits<G>::site_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = sites(g); vi != vi_end; ++vi) {
      if (site(*vi, g).s.get_twice() > 1 &&
          alps::is_positive<1>(site(*vi, g).d)) {
        frustrated = true;
        break;
      }
    }
  }
  return frustrated;
}

} // end namespace looper

#endif // LOOPER_MODEL_H
