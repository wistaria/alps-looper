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

#include <looper/lapack.h>
#include <looper/operator.h>
#include <looper/util.h>

#include <alps/parameterlist.h>
#include <alps/lattice.h>
#include <alps/model.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

template<class VEC>
std::pair<double, double> static_average2(double beta, double offset,
                                          const VEC& evals) {
  typedef VEC vector_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  double val = 0.0;
  double val2 = 0.0;
  for (; eval != eval_end; ++eval) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    val += (*eval) * weight;
    val2 += looper::sqr(*eval) * weight;
  }
  return std::make_pair(val, val2);
}

template<class VEC, class MAT>
double static_average(double beta, double offset,
                      const VEC& evals, const MAT& evecs,
                      const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  typename matrix_type::const_reverse_iterator1 evec = evecs.rbegin1();
  double val = 0.0;
  for (; eval != eval_end; ++eval, ++evec) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    typename matrix_type::const_iterator2 j = evec.begin();
    typename vector_type::const_iterator op = diagonal_matrix.begin();
    double v = 0.0;
    for (; j != evec.end(); ++j, ++op) v += looper::sqr(*j) * (*op);
    val += v * weight;
  }
  return val;
}

template<class VEC, class MAT>
std::pair<double, double> static_average2(double beta, double offset,
                                          const VEC& evals, const MAT& evecs,
                                          const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval = evals.rbegin();
  typename vector_type::const_reverse_iterator eval_end = evals.rend();
  typename matrix_type::const_reverse_iterator1 evec = evecs.rbegin1();
  double val = 0.0;
  double val2 = 0.0;
  for (; eval != eval_end; ++eval, ++evec) {
    double weight = std::exp(- beta * (*eval - offset)); // Boltzman weight
    typename matrix_type::const_iterator2 j = evec.begin();
    typename vector_type::const_iterator op = diagonal_matrix.begin();
    double v = 0.0;
    double v2 = 0.0;
    for (; j != evec.end(); ++j, ++op) {
      v += looper::sqr(*j) * (*op);
      v2 += looper::sqr(*j) * looper::sqr(*op);
    }
    val += v * weight;
    val2 += v2 * weight;
  }
  return std::make_pair(val, val2);
}

template<class VEC, class MAT>
double dynamic_average2(double beta, double offset,
                        const VEC& evals, const MAT& evecs,
                        const VEC& diagonal_matrix) {
  typedef VEC vector_type;
  typedef MAT matrix_type;
  typename vector_type::const_reverse_iterator eval0 = evals.rbegin();
  typename vector_type::const_reverse_iterator eval0_end = evals.rend();
  typename matrix_type::const_reverse_iterator1 evec0 = evecs.rbegin1();
  double val = 0.0;
  for (; eval0 != eval0_end; ++eval0, ++evec0) {
    double weight = std::exp(- beta * (*eval0 - offset)); // Boltzman weight
    typename vector_type::const_reverse_iterator eval1 = evals.rbegin();
    typename matrix_type::const_reverse_iterator1 evec1 = evecs.rbegin1();
    for (; evec1 != evec0; ++eval1, ++evec1) {
      // for evec0 != evec1
      double v = 0.;
      typename matrix_type::const_iterator2 v0 = evec0.begin();
      typename matrix_type::const_iterator2 v0_end = evec0.end();
      typename matrix_type::const_iterator2 v1 = evec1.begin();
      typename vector_type::const_iterator op = diagonal_matrix.begin();
      for (; v0 != v0_end; ++v0, ++v1, ++op) v += (*v0) * (*op) * (*v1);
      double wij;
      if (std::abs(*eval0 - *eval1) > 1.e-12) {
        wij = - (weight - exp(- beta * (*eval1 - offset))) / (*eval0 - *eval1);
      } else {
        wij = beta * weight;
      }
      val += 2 * looper::sqr(v) * wij;
    }
    {
      // for evec0 = evec1
      double v = 0.;
      typename matrix_type::const_iterator2 v0 = evec0.begin();
      typename matrix_type::const_iterator2 v0_end = evec0.end();
      typename matrix_type::const_iterator2 v1 = evec1.begin();
      typename vector_type::const_iterator op = diagonal_matrix.begin();
      for (; v0 != v0_end; ++v0, ++v1, ++op) v += (*v0) * (*op) * (*v1);
      double wij = beta * weight;
      val += looper::sqr(v) * wij;
    }
  }
  return val;
}


int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  typedef boost::numeric::ublas::vector<double> diagonal_matrix_type;

  using looper::sqr;

  alps::ParameterList parameterlist;
  std::cin >> parameterlist;
  for (alps::ParameterList::const_iterator p = parameterlist.begin();
       p != parameterlist.end(); ++p) {

    //
    // parameters
    //

    alps::Parameters params(*p);
    assert(params.defined("T"));
    for (alps::Parameters::const_iterator ps = p->begin(); ps != p->end(); ++ps)
      if (ps->key() != "LATTICE_LIBRARY" && ps->key() != "MODEL_LIBRARY")
        std::cout << ps->key() << " = " << ps->value() << std::endl;
    double beta = 1.0 / static_cast<double>(params["T"]);

    //
    // lattice & graph
    //

    typedef alps::graph_helper<>::graph_type graph_type;
    alps::graph_helper<> lattice(params);
    int num_vertices = boost::num_vertices(lattice.graph());
    bool is_bipartite = alps::set_parity(lattice.graph());

    //
    // model
    //

    alps::model_helper<> model(params);
    params.copy_undefined(model.model().default_parameters());

    //
    // generate basis set
    //

    alps::basis_states<short>
      basis_set(alps::basis_states_descriptor<short>(model.basis(),
                                                   lattice.graph()));
    int dim = basis_set.size();
    std::cout << "dimension of matrix = " << dim << std::endl;

    //
    // generate Hamiltonian matrix
    //

    matrix_type hamiltonian(dim, dim);
    hamiltonian.clear();
    alps::graph_traits<graph_type>::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(lattice.graph());
         vi != vi_end; ++vi)
      looper::add_to_matrix(hamiltonian,
        model.model().site_term(boost::get(alps::vertex_type_t(),
                                           lattice.graph(), *vi)),
        model.operators(), basis_set, *vi, lattice.graph(), params);
    alps::graph_traits<graph_type>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(lattice.graph());
         ei != ei_end; ++ei)
      looper::add_to_matrix(hamiltonian,
        model.model().bond_term(boost::get(alps::edge_type_t(),
                                           lattice.graph(), *ei)),
        model.operators(), basis_set,
        boost::source(*ei, lattice.graph()),
        boost::target(*ei, lattice.graph()),
        lattice.graph(), params);

    //
    // diagonalization
    //

    vector_type evals(dim);
    std::cout << "diagonalization... " << std::flush;
    looper::diagonalize(hamiltonian, evals);
    std::cout << "done\n";

    //
    // partition function, energy and specific heat
    //

    double gs_ene = evals(0);
    double part = 0.;
    vector_type::reverse_iterator eval_end = evals.rend();
    for (vector_type::reverse_iterator eval = evals.rbegin();
         eval != eval_end; ++eval) {
      double weight = std::exp(- beta * (*eval - gs_ene)); // Boltzman weight
      part += weight; // partition function
    }

    double ene, ene2;
    boost::tie(ene, ene2) = static_average2(beta, gs_ene, evals);
    ene = ene / part / num_vertices;
    ene2 = ene2 / part / sqr(num_vertices);
    double c = sqr(beta) * num_vertices * (ene2 - sqr(ene));

    std::cout << "ground state energy          = " << gs_ene << std::endl
              << "ground state energy per site = "
              << gs_ene/num_vertices << std::endl
              << "energy per site              = " << ene << std::endl
              << "specific heat                = " << c << std::endl;

    //
    // generate uniform/staggered Sz matrix
    //

    diagonal_matrix_type uniform_sz(dim);
    uniform_sz.clear();
    for (boost::tie(vi, vi_end) = boost::vertices(lattice.graph());
         vi != vi_end; ++vi) {
      looper::add_to_diagonal_matrix(uniform_sz,
        alps::SiteTermDescriptor<short>("Sz"), model.operators(),
        basis_set, *vi, lattice.graph());
    }

    double umag, umag2;
    boost::tie(umag, umag2) =
      static_average2(beta, gs_ene, evals, hamiltonian, uniform_sz);
    if (std::abs(umag) < 1.0e-12) umag = 0.;
    double usus = dynamic_average2(beta, gs_ene, evals, hamiltonian, uniform_sz);
    umag = umag / part / num_vertices;
    umag2 = umag2 / part / num_vertices;
    usus = usus / part / num_vertices;
    std::cout << "uniform magnetization        = " << umag << std::endl
              << "uniform magnetization^2      = " << umag2 << std::endl
              << "uniform susceptibility       = " << usus << std::endl;

    if (is_bipartite) {
      diagonal_matrix_type staggered_sz(dim);
      staggered_sz.clear();
      for (boost::tie(vi, vi_end) = boost::vertices(lattice.graph());
           vi != vi_end; ++vi) {
        if (boost::get(alps::parity_t(), lattice.graph(), *vi) ==
            alps::parity::white) {
          looper::add_to_diagonal_matrix(staggered_sz,
            alps::SiteTermDescriptor<short>("Sz"), model.operators(),
            basis_set, *vi, lattice.graph());
        } else {
          looper::add_to_diagonal_matrix(staggered_sz,
            alps::SiteTermDescriptor<short>("-Sz"), model.operators(),
            basis_set, *vi, lattice.graph());
        }
      }

      double smag, smag2;
      boost::tie(smag, smag2) =
        static_average2(beta, gs_ene, evals, hamiltonian, staggered_sz);
      if (std::abs(smag) < 1.0e-12) smag = 0.;
      double ssus = dynamic_average2(beta, gs_ene, evals, hamiltonian, staggered_sz);
      smag = smag / part / num_vertices;
      smag2 = smag2 / part / num_vertices;
      ssus = ssus / part / num_vertices;
      std::cout << "staggered magnetization      = " << smag << std::endl
                << "staggered magnetization^2    = " << smag2 << std::endl
                << "staggered susceptibility     = " << ssus << std::endl;
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}
