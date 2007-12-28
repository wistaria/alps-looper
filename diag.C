/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2003-2007 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_config.h"
#include "loop_factory.h"
#include <looper/lapack.h>
#include <looper/power.h>
#include <alps/math.hpp>
#include <alps/scheduler.h>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/tuple/tuple.hpp>

namespace {

template<class MATRIX, class I, class GRAPH>
void add_to_matrix(
  MATRIX& matrix,
  const alps::HamiltonianDescriptor<I>& hd,
  const alps::BasisDescriptor<I>& basis,
  const alps::basis_states<I>& basis_set,
  const typename alps::graph_traits<GRAPH>::vertex_descriptor vd,
  const GRAPH& graph,
  const alps::Parameters& params /* = alps::Parameters() */) {

  typedef typename MATRIX::value_type value_type;
  typedef alps::basis_states<I> basis_set_type;

  int t = get(alps::site_type_t(), graph, vd);
  int s = get(alps::site_index_t(), graph, vd);
  int dim = basis_set.size();
  int ds = basis_set.basis().get_site_basis(s).num_states();

  boost::multi_array<value_type, 2>
    site_matrix(get_matrix(value_type(), hd.site_term(t), basis.site_basis(t), params));

  for (int i = 0; i < dim; ++i) {
    int is = basis_set[i][s];
    for (int js = 0; js < ds; ++js) {
      typename basis_set_type::value_type target(basis_set[i]);
      target[s] = js;
      int j = basis_set.index(target);
      if (j < dim) matrix(i,j) += site_matrix[is][js];
    }
  }
}

template<class MATRIX, class I, class GRAPH>
void add_to_matrix(
  MATRIX& matrix,
  const alps::HamiltonianDescriptor<I>& hd,
  const alps::BasisDescriptor<I>& basis,
  const alps::basis_states<I>& basis_set,
  const typename alps::graph_traits<GRAPH>::bond_descriptor ed,
  const typename alps::graph_traits<GRAPH>::site_descriptor vd0,
  const typename alps::graph_traits<GRAPH>::site_descriptor vd1,
  const GRAPH& graph,
  const alps::Parameters& params /* = alps::Parameters() */) {

  typedef typename MATRIX::value_type value_type;
  typedef alps::basis_states<I> basis_set_type;

  int t = get(alps::bond_type_t(), graph, ed);
  int st0 = get(alps::site_type_t(), graph, vd0);
  int st1 = get(alps::site_type_t(), graph, vd1);
  int s0 = get(alps::site_index_t(), graph, vd0);
  int s1 = get(alps::site_index_t(), graph, vd1);
  int dim = basis_set.size();
  int ds0 = basis_set.basis().get_site_basis(s0).num_states();
  int ds1 = basis_set.basis().get_site_basis(s1).num_states();

  boost::multi_array<value_type, 4>
    bond_matrix(alps::get_matrix(
      value_type(), hd.bond_term(t),
      basis.site_basis(st0), basis.site_basis(st1),
      params));

  for (int i = 0; i < dim; ++i) {
    int is0 = basis_set[i][s0];
    int is1 = basis_set[i][s1];
    for (int js0 = 0; js0 < ds0; ++js0) {
      for (int js1 = 0; js1 < ds1; ++js1) {
        typename basis_set_type::value_type target(basis_set[i]);
        target[s0] = js0;
        target[s1] = js1;
        int j = basis_set.index(target);
        if (j < dim) matrix(i,j) += bond_matrix[is0][is1][js0][js1];
      }
    }
  }
}

template<class VECTOR, class I, class GRAPH>
void add_to_diagonal_matrix(
  VECTOR& vector,
  const alps::SiteTermDescriptor& term,
  const alps::BasisDescriptor<I>& basis,
  const alps::basis_states<I>& basis_set,
  const typename alps::graph_traits<GRAPH>::vertex_descriptor& vd,
  const GRAPH& graph,
  const alps::Parameters& params /* = alps::Parameters() */) {

  typedef typename VECTOR::value_type value_type;

  int t = get(alps::site_type_t(), graph, vd);
  int s = get(alps::site_index_t(), graph, vd);
  int dim = basis_set.size();
  int ds = basis_set.basis().get_site_basis(s).num_states();

  boost::multi_array<value_type, 2>
    site_matrix(get_matrix(value_type(), term, basis.site_basis(t), params));
  for (int is = 0; is < ds; ++is)
    for (int js = 0; js < ds; ++js)
      if ((is != js) && alps::is_nonzero<1>(site_matrix[is][js]))
        boost::throw_exception(std::logic_error("non-diagonal site term"));

  for (int i = 0; i < dim; ++i) {
    int is = basis_set[i][s];
    vector(i) += site_matrix[is][is];
  }
}

template<class VEC>
std::pair<double, double> static_average2(double beta, double offset, const VEC& evals) {
  typedef VEC vector_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  double val = 0.0;
  double val2 = 0.0;
  for (; eval != eval_end; ++eval) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    val += (*eval) * weight;
    val2 += looper::power2(*eval) * weight;
  }
  return std::make_pair(val, val2);
}

template<class VEC, class MAT>
double static_average(double beta, double offset, const VEC& evals, const MAT& evecs,
                      const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  typename matrix_type::const_reverse_iterator2 evec = evecs.rbegin2();
  double val = 0.0;
  for (; eval != eval_end; ++eval, ++evec) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    typename matrix_type::const_iterator1 j = evec.begin();
    typename vector_type::const_iterator op = diagonal_matrix.begin();
    double v = 0.0;
    for (; j != evec.end(); ++j, ++op) v += looper::power2(*j) * (*op);
    val += v * weight;
  }
  return val;
}

template<class VEC, class MAT>
std::pair<double, double> static_average2(double beta, double offset, const VEC& evals,
                                          const MAT& evecs, const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  typename matrix_type::const_reverse_iterator2 evec = evecs.rbegin2();
  double val = 0.0;
  double val2 = 0.0;
  for (; eval != eval_end; ++eval, ++evec) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    typename matrix_type::const_iterator1 j = evec.begin();
    typename vector_type::const_iterator op = diagonal_matrix.begin();
    double v = 0.0;
    double v2 = 0.0;
    for (; j != evec.end(); ++j, ++op) {
      v += looper::power2(*j) * (*op);
      v2 += looper::power2(*j) * looper::power2(*op);
    }
    val += v * weight;
    val2 += v2 * weight;
  }
  return std::make_pair(val, val2);
}

template<class VEC, class MAT>
boost::tuple<double, double, double>
static_average4(double beta, double offset, const VEC& evals, const MAT& evecs,
                const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  typename matrix_type::const_reverse_iterator2 evec = evecs.rbegin2();
  double val = 0.0;
  double val2 = 0.0;
  double val4 = 0.0;
  for (; eval != eval_end; ++eval, ++evec) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    typename matrix_type::const_iterator1 j = evec.begin();
    typename vector_type::const_iterator op = diagonal_matrix.begin();
    double v = 0.0;
    double v2 = 0.0;
    double v4 = 0.0;
    for (; j != evec.end(); ++j, ++op) {
      v += looper::power2(*j) * (*op);
      v2 += looper::power2(*j) * looper::power2(*op);
      v4 += looper::power2(*j) * looper::power4(*op);
    }
    val += v * weight;
    val2 += v2 * weight;
    val4 += v4 * weight;
  }
  return boost::make_tuple(val, val2, val4);
}

template<class VEC, class MAT>
double dynamic_average2(double beta, double offset, const VEC& evals, const MAT& evecs,
                        const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval0 = evals.rbegin();
  typename vector_type::const_reverse_iterator eval0_end = evals.rend();
  typename matrix_type::const_reverse_iterator2 evec0 = evecs.rbegin2();
  double val = 0.0;
  for (; eval0 != eval0_end; ++eval0, ++evec0) {
    double weight = std::exp(- beta * (*eval0 - offset)); // Boltzman weight
    typename vector_type::const_reverse_iterator eval1 = evals.rbegin();
    typename matrix_type::const_reverse_iterator2 evec1 = evecs.rbegin2();
    for (; evec1 != evec0; ++eval1, ++evec1) {
      // for evec0 != evec1
      double v = 0.;
      typename matrix_type::const_iterator1 v0 = evec0.begin();
      typename matrix_type::const_iterator1 v0_end = evec0.end();
      typename matrix_type::const_iterator1 v1 = evec1.begin();
      typename vector_type::const_iterator op = diagonal_matrix.begin();
      for (; v0 != v0_end; ++v0, ++v1, ++op) v += (*v0) * (*op) * (*v1);
      double wij;
      if (std::abs(*eval0 - *eval1) > 1.e-12) {
        wij = - (weight - exp(- beta * (*eval1 - offset))) / (*eval0 - *eval1);
      } else {
        wij = beta * weight;
      }
      val += 2 * looper::power2(v) * wij;
    }
    {
      // for evec0 = evec1
      double v = 0.;
      typename matrix_type::const_iterator1 v0 = evec0.begin();
      typename matrix_type::const_iterator1 v0_end = evec0.end();
      typename matrix_type::const_iterator1 v1 = evec1.begin();
      typename vector_type::const_iterator op = diagonal_matrix.begin();
      for (; v0 != v0_end; ++v0, ++v1, ++op) v += (*v0) * (*op) * (*v1);
      double wij = beta * weight;
      val += looper::power2(v) * wij;
    }
  }
  return val;
}


class diag_worker :
  protected looper::lattice_helper<loop_config::lattice_graph_t>,
  protected alps::model_helper<> {
public:
  diag_worker(alps::Parameters const& p, alps::ObservableSet&) :
    looper::lattice_helper<loop_config::lattice_graph_t>(p),
    alps::model_helper<>(this->graph_helper(), p), done(false), params(p) {}
  bool is_thermalized() const { return true; }
  double progress() const { return done ? 1 : 0; }
  template<typename ENGINE>
  void run(ENGINE&, alps::ObservableSet& obs);
  void set_beta(double) const { boost::throw_exception(std::logic_error("diag")); }
  double g_weight() const { boost::throw_exception(std::logic_error("diag")); return 0; }
  double lambda(double) const { boost::throw_exception(std::logic_error("diag")); return 0; }
  void save(alps::ODump& dp) const { dp << done; }
  void load(alps::IDump& dp) { dp >> done; }
private:
  bool done;
  alps::Parameters params;
};


template<typename ENGINE>
void diag_worker::run(ENGINE&, alps::ObservableSet& obs) {
  typedef typename looper::real_site_descriptor<
    typename looper::lattice_helper<loop_config::lattice_graph_t> >::type site_descriptor;
  typedef typename looper::real_bond_descriptor<
    typename looper::lattice_helper<loop_config::lattice_graph_t> >::type bond_descriptor;

  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
  typedef boost::numeric::ublas::vector<double> diagonal_matrix_type;

  if (done) return;
  done = true;

  // parameters
  double beta = 1.0 / alps::evaluate("T", params);
  if (beta < 0)
    boost::throw_exception(std::invalid_argument("negative temperature"));
  double nsite = this->graph_helper().num_sites();
  double volume = this->volume();

  // measurements
  std::map<std::string, double> m;
  m["Temperature"] = 1/beta;
  m["Inverse Temperature"] = beta;
  m["Number of Sites"] = nsite;
  m["Volume"] = volume;

  // generate basis set
  alps::basis_states<short>
    basis_set(alps::basis_states_descriptor<short>(model().basis(), this->graph_helper().graph()));
  int dim = basis_set.size();
  m["Dimension of Matrix"] = dim;

  // generate Hamiltonian matrix
  matrix_type hamiltonian(dim, dim);
  hamiltonian.clear();
  BOOST_FOREACH(site_descriptor v, this->graph_helper().sites())
    add_to_matrix(hamiltonian, model(), model().basis(), basis_set, v, this->rg(), params);
  BOOST_FOREACH(bond_descriptor e, this->graph_helper().bonds())
    add_to_matrix(hamiltonian, model(), model().basis(), basis_set, e, source(e, this->rg()),
    target(e, this->rg()), this->rg(), params);
  diagonal_matrix_type diagonal_energy(dim);
  for (int i = 0; i < dim; ++i) diagonal_energy(i) = hamiltonian(i,i);

  // partition function coefficients
  matrix_type mat(dim, dim);
  mat.clear();
  for (int i = 0; i < dim; ++i) mat(i,i) = 1;
  double f = 1;
  for (int k = 0; k < 5; ++k) {
    if (k != 0) {
      mat = prod(hamiltonian, mat);
      f *= -k;
    }
    double tr = 0;
    for (int i = 0; i < dim; ++i) tr += mat(i,i);
    m["Partition Function Coefficient #" + boost::lexical_cast<std::string>(k)] = tr / f;
  }

  // diagonalization
  vector_type evals(dim);
  std::cerr << "start diagonalization... " << std::flush;
  looper::diagonalize(hamiltonian, evals);
  std::cerr << "done\n";

  // partition function, energy and specific heat
  double gs_ene = evals(0);
  double part = 0.;
  vector_type::reverse_iterator eval_end = evals.rend();
  for (vector_type::reverse_iterator eval = evals.rbegin(); eval != eval_end; ++eval) {
    double weight = std::exp(- beta * (*eval - gs_ene)); // Boltzmann weight
    part += weight; // partition function
  }

  double ene, ene2;
  boost::tie(ene, ene2) = static_average2(beta, gs_ene, evals);
  ene = ene / part;
  ene2 = ene2 / part / looper::power2(volume);
  double fe = gs_ene - std::log(part) / beta;
  m["Ground State Energy"] = gs_ene;
  m["Ground State Energy Density"] = gs_ene / volume;
  m["Energy"] = ene;
  m["Energy Density"] = ene / volume;
  m["Specific Heat"] = looper::power2(beta) * volume * (ene2 - looper::power2(ene/volume));
  m["Free Energy"] = fe;
  m["Free Energy Density"] = fe / volume;
  m["Entropy"] = beta * (ene - fe);
  m["Entropy Density"] = beta * (ene - fe) / volume;

  // generate uniform/staggered Sz matrix
  diagonal_matrix_type uniform_sz(dim);
  uniform_sz.clear();
  BOOST_FOREACH(site_descriptor v, this->graph_helper().sites())
    add_to_diagonal_matrix(uniform_sz, alps::SiteTermDescriptor("Sz(i)", "i"),
                           model().basis(), basis_set, v, this->rg(), params);
  double umag, umag2, umag4;
  boost::tie(umag, umag2, umag4) = static_average4(beta, gs_ene, evals, hamiltonian, uniform_sz);
  umag = alps::round<1>(umag);
  m["Magnetization"] = umag / part;
  m["Magnetization Density"] = umag / part / volume;
  m["Magnetization^2"] = umag2 / part;
  m["Magnetization Density^2"] = umag2 / part / looper::power2(volume);
  m["Magnetization^4"] = umag4 / part;
  m["Magnetization Density^4"] = umag4 / part / looper::power4(volume);
  m["Susceptibility"] =
    dynamic_average2(beta, gs_ene, evals, hamiltonian, uniform_sz) / part / volume;
  m["Binder Ratio of Magnetization"] = looper::power2(umag2 / part) / (umag4 / part);
  if (this->graph_helper().is_bipartite()) {
    diagonal_matrix_type staggered_sz(dim);
    staggered_sz.clear();
    BOOST_FOREACH(site_descriptor v, this->graph_helper().sites()) {
      int g = 2 * get(alps::parity_t(), this->rg(), v) - 1;
      if (g == 1) {
        add_to_diagonal_matrix(staggered_sz, alps::SiteTermDescriptor("Sz(i)", "i"),
                               model().basis(), basis_set, v, this->rg(), params);
      } else {
        add_to_diagonal_matrix(staggered_sz, alps::SiteTermDescriptor("-Sz(i)", "i"),
                               model().basis(), basis_set, v, this->rg(), params);
      }
    }
    double smag, smag2, smag4;
    boost::tie(smag, smag2, smag4) =
      static_average4(beta, gs_ene, evals, hamiltonian, staggered_sz);
    smag = alps::round<1>(smag);
    m["Staggered Magnetization"] = smag / part;
    m["Staggered Magnetization Density"] = smag / part / volume;
    m["Staggered Magnetization^2"] = smag2 / part;
    m["Staggered Magnetization Density^2"] = smag2 / part / looper::power2(volume);
    m["Staggered Magnetization^4"] = smag4 / part;
    m["Staggered Magnetization Density^4"] = smag4 / part / looper::power4(volume);
    m["Staggered Susceptibility"] =
      dynamic_average2(beta, gs_ene, evals, hamiltonian, staggered_sz) / part / volume;
    m["Binder Ratio of Staggered Magnetization"] = looper::power2(smag2 / part) / (smag4 / part);

    // local Sz
    diagonal_matrix_type local_sz(dim);
    BOOST_FOREACH(site_descriptor v, this->graph_helper().sites()) {
      local_sz.clear();
      add_to_diagonal_matrix(local_sz, alps::SiteTermDescriptor("Sz(i)", "i"),
                             model().basis(), basis_set, v, this->rg(), params);
      double mag, mag2, mag4;
      boost::tie(mag, mag2, mag4) = static_average4(beta, gs_ene, evals, hamiltonian, local_sz);
      mag = alps::round<1>(mag);
      m["Local Magnetization at Site " + boost::lexical_cast<std::string>(v)] = mag / part;
    }
  }

  // store measurements
  typedef std::pair<std::string, double> nv_t;
  BOOST_FOREACH(nv_t const& t, m) obs << alps::SimpleRealObservable(t.first);
  obs.reset(true);
  BOOST_FOREACH(nv_t const& t, m) obs[t.first] << t.second;
}


class dummy_evaluator : public looper::abstract_evaluator {
public:
  void pre_evaluate(alps::ObservableSet&, alps::Parameters const&,
    alps::ObservableSet const&) const {}
  void evaluate(alps::scheduler::MCSimulation&, alps::Parameters const&,
    boost::filesystem::path const&) const {}
  void evaluate(alps::ObservableSet&, alps::Parameters const&, alps::ObservableSet const&) const {}
};


//
// dynamic registration to the factories
//

const bool worker_registered = loop_factory::instance()->
  register_worker<diag_worker>("diagonalization");
const bool evaluator_registered = loop_factory::instance()->
  register_evaluator<dummy_evaluator>("diagonalization");

} // end namespace
