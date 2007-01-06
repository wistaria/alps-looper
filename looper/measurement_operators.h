/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_MEASUREMENT_OPERATORS_H
#define LOOPER_MEASUREMENT_OPERATORS_H

#include <alps/model.h>
#include <alps/scheduler/measurement_operators.h>
#include <boost/foreach.hpp>
#include <boost/next_prior.hpp>

namespace looper {

class measurement_operators : private alps::MeasurementOperators {
public:
  typedef alps::MeasurementOperators super_type;

  template<typename RG, typename I>
  measurement_operators(alps::Parameters const& params, lattice_helper<RG> const& lat,
    alps::model_helper<I> const& mh) : alps::MeasurementOperators(params) {
    using std::max;
    max_type_ = 0;
    BOOST_FOREACH(typename real_site_descriptor<lattice_helper<RG> >::type const& s,
      sites(lat.rg())) {
      unsigned int t = get(site_type_t(), lat.rg(), s);
      types_.insert(t);
      max_type_ = max(max_type_, t);
    }

    for (std::map<std::string, std::string>::iterator itr =
         super_type::average_expressions.begin();
         itr != super_type::average_expressions.end(); /* no increment */) {
      ++itr;
      if (!build_diagonal_matrix(boost::prior(itr)->second, params, mh)) {
        std::clog << "Will not measure \"" << boost::prior(itr)->first
                  << "\" since it is off-diagonal\n";
        super_type::average_expressions.erase(boost::prior(itr));
      }
    }
    for (std::map<std::string, std::string>::iterator itr = super_type::local_expressions.begin();
         itr != super_type::local_expressions.end(); /* no increment */) {
      ++itr;
      if (!build_diagonal_matrix(boost::prior(itr)->second, params, mh)) {
        std::clog << "Will not measure \"" << boost::prior(itr)->first
                  << "\" since it is off-diagonal\n";
        super_type::local_expressions.erase(boost::prior(itr));
      }
    }
    for (std::map<std::string, std::pair<std::string, std::string> >::iterator
         itr = super_type::correlation_expressions.begin();
         itr != super_type::correlation_expressions.end(); /* no increment */) {
      ++itr;
      if (!build_diagonal_matrix(boost::prior(itr)->second.first, params, mh) ||
          !build_diagonal_matrix(boost::prior(itr)->second.second, params, mh)) {
        std::clog << "Will not measure \"" << boost::prior(itr)->first
                  << "\" since it is off-diagonal\n";
        super_type::correlation_expressions.erase(boost::prior(itr));
      }
    }
    for (std::map<std::string, std::pair<std::string, std::string> >::iterator
         itr = super_type::structurefactor_expressions.begin();
         itr != super_type::structurefactor_expressions.end(); /* no increment */) {
      ++itr;
      if (!build_diagonal_matrix(boost::prior(itr)->second.first, params, mh) ||
          !build_diagonal_matrix(boost::prior(itr)->second.second, params, mh)) {
        std::clog << "Will not measure \"" << boost::prior(itr)->first
                  << "\" since it is off-diagonal\n";
        super_type::structurefactor_expressions.erase(boost::prior(itr));
      }
    }
  }

protected:
  template<class I>
  bool build_diagonal_matrix(std::string const& op, alps::Parameters const& p,
    alps::model_helper<I> const& mh) {
    alps::SiteOperator term(op);
    bool valid = true;
    std::vector<std::vector<double> > mat(max_type_+1);
    BOOST_FOREACH(unsigned int t, types_) {
      boost::multi_array<double, 2> m =
        get_matrix(double(), term, mh.model().basis().site_basis(t), p);
      std::size_t dim = m.shape()[0];
      mat[t].resize(dim);
      for (std::size_t i = 0; i < dim; ++i)
        for (std::size_t j = 0; j < dim; ++j)
          if (i == j)
            mat[t][i] = m[i][j];
          else
            if (alps::is_nonzero(m[i][j])) valid = false;
    }
    if (valid) diagonal_matrix_element_[op] = mat;
    return valid;
  }

private:
  std::set<unsigned int> types_;
  unsigned int max_type_;
  std::map<std::string, std::vector<std::vector<double> > > diagonal_matrix_element_;
};

} // end namespace looper

#endif // LOOPER_MEASUREMENT_OPERATORS_H
