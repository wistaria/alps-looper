/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2003-2004 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_OPERATOR_H
#define LOOPER_OPERATOR_H

#include <alps/model.h>

namespace looper {

typedef alps::SiteBasisDescriptor<short> site_basis_descriptor;
typedef std::map<std::string, alps::OperatorDescriptor<short> >
  operator_map_type;

template<class I>
site_basis_descriptor spin_basis(const alps::half_integer<I>& s = 0.5)
{
  alps::Parameters p;
  p["S"] = s;
  site_basis_descriptor sb("spin", p);
  sb.push_back(alps::QuantumNumberDescriptor<short>("Sz", "-S", "S"));
  sb.set_parameters(alps::Parameters()); // required for re-evaluation
  return sb;
}

inline site_basis_descriptor spin_basis(double s)
{ return spin_basis(alps::half_integer<short>(s)); }

inline operator_map_type spin_operators()
{
  operator_map_type ops;

  // identity operator
  ops["I"] = alps::OperatorDescriptor<short>("I", "1");

  // Sz
  ops["Sz"] = alps::OperatorDescriptor<short>("Sz", "Sz");

  // S+
  alps::OperatorDescriptor<short> splus("Splus", "sqrt(S*(S+1)-Sz*(Sz+1))");
  splus["Sz"] = 1;
  ops["Splus"] = splus;

  // S-
  alps::OperatorDescriptor<short> sminus("Sminus", "sqrt(S*(S+1)-Sz*(Sz-1))");
  sminus["Sz"] = -1;
  ops["Sminus"] = sminus;

  return ops;
}

template<class MATRIX, class I, class GRAPH>
void add_to_matrix(MATRIX& matrix,
                   const alps::SiteTermDescriptor<I>& term,
                   const operator_map_type& ops,
                   const alps::basis_states<I>& basis_states,
                   const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd,
                   const GRAPH& graph,
                   const alps::Parameters& params = alps::Parameters())
{
  typedef alps::basis_states<I> basis_state_type;

  int s = boost::get(boost::vertex_index_t(), graph, vd);
  int dim = basis_states.size();
  int ds = basis_states.basis().get_site_basis(s).num_states();

  boost::multi_array<double, 2>
    site_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s),
                                 ops,
                                   params));

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is = state[s];
    for (int js = 0; js < ds; ++js) {
      if (site_matrix[is][js] != 0.0) {
        typename basis_state_type::value_type target = state;
        target[s] = js;
        int j = basis_states.index(target);
        if (j < dim) matrix(i,j) += site_matrix[is][js];
      }
    }
  }
}

template<class VECTOR, class I, class GRAPH>
void apply_to_vector(const VECTOR& vec_in, VECTOR& vec_out,
                     const alps::SiteTermDescriptor<I>& term,
                     const operator_map_type& ops,
                     const alps::basis_states<I>& basis_states,
                     const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd,
                     const GRAPH& graph,
                     const alps::Parameters& params = alps::Parameters())
{
  typedef alps::basis_states<I> basis_state_type;

  int s = boost::get(boost::vertex_index_t(), graph, vd);
  int dim = basis_states.size();
  int ds = basis_states.basis().get_site_basis(s).num_states();

  boost::multi_array<double, 2>
    site_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s),
                                 ops,
                                   params));

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is = state[s];
    for (int js = 0; js < ds; ++js) {
      if (site_matrix[is][js] != 0.0) {
        typename basis_state_type::value_type target = state;
        target[s] = js;
        int j = basis_states.index(target);
        if (j < dim) vec_out(j) += site_matrix[is][js] * vec_in(i);
      }
    }
  }
}

template<class MATRIX, class I, class GRAPH>
void add_to_matrix(MATRIX& matrix,
                   const alps::BondTermDescriptor<I>& term,
                   const operator_map_type& ops,
                   const alps::basis_states<I>& basis_states,
                   const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd0,
                   const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd1,
                   const GRAPH& graph,
                   const alps::Parameters& params = alps::Parameters())
{
  typedef alps::basis_states<I> basis_state_type;

  int s0 = boost::get(boost::vertex_index_t(), graph, vd0);
  int s1 = boost::get(boost::vertex_index_t(), graph, vd1);
  int dim = basis_states.size();
  int ds0 = basis_states.basis().get_site_basis(s0).num_states();
  int ds1 = basis_states.basis().get_site_basis(s1).num_states();

  boost::multi_array<double, 4>
    bond_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s0),
                                 basis_states.basis().get_site_basis(s1),
                                 ops,
                                 params));

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is0 = state[s0];
    int is1 = state[s1];
    for (int js0 = 0; js0 < ds0; ++js0) {
      for (int js1 = 0; js1 < ds1; ++js1) {
        if (bond_matrix[is0][is1][js0][js1] != 0.0) {
          typename basis_state_type::value_type target = state;
          target[s0] = js0;
          target[s1] = js1;
          int j = basis_states.index(target);
          if (j < dim) matrix(i,j) += bond_matrix[is0][is1][js0][js1];
        }
      }
    }
  }
}

template<class VECTOR, class I, class GRAPH>
void apply_to_vector(const VECTOR& vec_in, VECTOR& vec_out,
                     const alps::BondTermDescriptor<I>& term,
                     const operator_map_type& ops,
                     const alps::basis_states<I>& basis_states,
                     const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd0,
                     const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd1,
                     const GRAPH& graph,
                     const alps::Parameters& params = alps::Parameters())
{
  typedef alps::basis_states<I> basis_state_type;

  int s0 = boost::get(boost::vertex_index_t(), graph, vd0);
  int s1 = boost::get(boost::vertex_index_t(), graph, vd1);
  int dim = basis_states.size();
  int ds0 = basis_states.basis().get_site_basis(s0).num_states();
  int ds1 = basis_states.basis().get_site_basis(s1).num_states();

  boost::multi_array<double, 4>
    bond_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s0),
                                 basis_states.basis().get_site_basis(s1),
                                 ops,
                                   params));

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is0 = state[s0];
    int is1 = state[s1];
    for (int js0 = 0; js0 < ds0; ++js0) {
      for (int js1 = 0; js1 < ds1; ++js1) {
        if (bond_matrix[is0][is1][js0][js1] != 0.0) {
          typename basis_state_type::value_type target = state;
          target[s0] = js0;
          target[s1] = js1;
          int j = basis_states.index(target);
          if (j < dim) vec_out(j) += bond_matrix[is0][is1][js0][js1] * vec_in(i);
        }
      }
    }
  }
}

template<class VECTOR, class I, class GRAPH>
void add_to_diagonal_matrix(VECTOR& vector,
                            const alps::SiteTermDescriptor<I>& term,
                            const operator_map_type& ops,
                            const alps::basis_states<I>& basis_states,
                            const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd,
                            const GRAPH& graph,
                            const alps::Parameters& params = alps::Parameters(),
                            double tol = 1.0e-10)
{
  typedef alps::basis_states<I> basis_state_type;

  int s = boost::get(boost::vertex_index_t(), graph, vd);
  int dim = basis_states.size();
  int ds = basis_states.basis().get_site_basis(s).num_states();

  boost::multi_array<double, 2>
    site_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s),
                                 ops,
                                   params));
  for (int is = 0; is < ds; ++is) {
    for (int js = 0; js < ds; ++js) {
      if ((is != js) && (std::abs(site_matrix[is][js]) > tol)) {
        boost::throw_exception(std::logic_error("non-diagonal site term"));
      }
    }
  }

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is = state[s];
    vector(i) += site_matrix[is][is];
  }
}

template<class VECTOR, class I, class GRAPH>
void apply_diagonal_to_vector(const VECTOR& vec_in, VECTOR& vec_out,
                              const alps::SiteTermDescriptor<I>& term,
                              const operator_map_type& ops,
                              const alps::basis_states<I>& basis_states,
                              const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd,
                              const GRAPH& graph,
                              const alps::Parameters& params = alps::Parameters(),
                              double tol = 1.0e-10)
{
  typedef alps::basis_states<I> basis_state_type;

  int s = boost::get(boost::vertex_index_t(), graph, vd);
  int dim = basis_states.size();
  int ds = basis_states.basis().get_site_basis(s).num_states();

  boost::multi_array<double, 2>
    site_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s),
                                 ops,
                                   params));
  for (int is = 0; is < ds; ++is) {
    for (int js = 0; js < ds; ++js) {
      if ((is != js) && (std::abs(site_matrix[is][js]) > tol)) {
        boost::throw_exception(std::logic_error("non-diagonal site term"));
      }
    }
  }

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is = state[s];
    vec_out(i) = site_matrix[is][is] * vec_in(i);
  }
}

template<class VECTOR, class I, class GRAPH>
void add_to_diagonal_matrix(VECTOR& vector,
                            const alps::BondTermDescriptor<I>& term,
                            const operator_map_type& ops,
                            const alps::basis_states<I>& basis_states,
                            const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd0,
                            const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd1,
                            const GRAPH& graph,
                            const alps::Parameters& params = alps::Parameters(),
                            double tol = 1.0e-10)
{
  typedef alps::basis_states<I> basis_state_type;

  int s0 = boost::get(boost::vertex_index_t(), graph, vd0);
  int s1 = boost::get(boost::vertex_index_t(), graph, vd1);
  int dim = basis_states.size();
  int ds0 = basis_states.basis().get_site_basis(s0).num_states();
  int ds1 = basis_states.basis().get_site_basis(s1).num_states();

  boost::multi_array<double, 4>
    bond_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s0),
                                 basis_states.basis().get_site_basis(s1),
                                 ops,
                                   params));
  for (int is0 = 0; is0 < ds0; ++is0) {
    for (int is1 = 0; is1 < ds1; ++is1) {
      for (int js0 = 0; js0 < ds0; ++js0) {
        for (int js1 = 0; js1 < ds1; ++js1) {
          if ((is0 != js0) && (is1 != js1) &&
              (std::abs(bond_matrix[is0][is1][js0][js1]) > tol)) {
            boost::throw_exception(std::logic_error("non-diagonal bond term"));
          }
        }
      }
    }
  }

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is0 = state[s0];
    int is1 = state[s1];
    vector(i) += bond_matrix[is0][is1][is0][is1];
  }
}

template<class VECTOR, class I, class GRAPH>
void apply_diagonal_to_vector(const VECTOR& vec_in, VECTOR& vec_out,
                              const alps::BondTermDescriptor<I>& term,
                              const operator_map_type& ops,
                              const alps::basis_states<I>& basis_states,
                              const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd0,
                              const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd1,
                              const GRAPH& graph,
                              const alps::Parameters& params = alps::Parameters(),
                              double tol = 1.0e-10)
{
  typedef alps::basis_states<I> basis_state_type;

  int s0 = boost::get(boost::vertex_index_t(), graph, vd0);
  int s1 = boost::get(boost::vertex_index_t(), graph, vd1);
  int dim = basis_states.size();
  int ds0 = basis_states.basis().get_site_basis(s0).num_states();
  int ds1 = basis_states.basis().get_site_basis(s1).num_states();

  boost::multi_array<double, 4>
    bond_matrix(alps::get_matrix(double(),
                                 term,
                                 basis_states.basis().get_site_basis(s0),
                                 basis_states.basis().get_site_basis(s1),
                                 ops,
                                   params));
  for (int is0 = 0; is0 < ds0; ++is0) {
    for (int is1 = 0; is1 < ds1; ++is1) {
      for (int js0 = 0; js0 < ds0; ++js0) {
        for (int js1 = 0; js1 < ds1; ++js1) {
          if ((is0 != js0) && (is1 != js1) &&
              (std::abs(bond_matrix[is0][is1][js0][js1]) > tol)) {
            boost::throw_exception(std::logic_error("non-diagonal bond term"));
          }
        }
      }
    }
  }

  for (int i = 0; i < dim; ++i) {
    typename basis_state_type::value_type state = basis_states[i];
    int is0 = state[s0];
    int is1 = state[s1];
    vec_out(i) += bond_matrix[is0][is1][is0][is1] * vec_in(i);
  }
}

} // end namespace looper

#endif // LOOPER_OPERATOR_H
