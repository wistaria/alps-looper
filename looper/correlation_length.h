/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_CORRELATION_LENGTH_H
#define LOOPER_CORRELATION_LENGTH_H

#include "measurement.h"
#include <boost/spirit/core.hpp>
#include <boost/spirit/actor.hpp>

namespace {

inline std::vector<double> parse_vec(std::string const& str, alps::Parameters const& params) {
  using namespace boost::spirit;
  std::vector<std::string> vec_str;
  subrule<0> vec;
  subrule<1> elem;
  if (!parse(str.c_str(),
             (vec  = '(' >> elem >> *(',' >> elem) >> ')',
              elem = (+(anychar_p - ',' - ')'))[push_back_a(vec_str)]
              ),
             space_p).full) {
    std::cerr << "can not parse: " << str << std::endl;
    boost::throw_exception(std::invalid_argument("correlation_length::parse"));
  }
  std::vector<double> result;
  for (int i = 0; i < vec_str.size(); ++i)
    result.push_back(2 * M_PI * alps::evaluate(vec_str[i], params));
  return result;
}

}

namespace looper {

struct correlation_length : public has_evaluator_tag {

  typedef dumb_measurement<correlation_length> dumb;

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef LAT lattice_t;
    typedef std::complex<double> complex_type;

    bool measure, improved;
    std::vector<complex_type> phase0, phase1;
    double dq_abs;
    std::vector<complex_type> sq0, sq1;

    template<typename M>
    void initialize(M& m, alps::Parameters const& params, lattice_t const& lat,
      bool is_signed, bool use_improved_estimator) {
      measure = params.value_or_default("MEASURE[Correlation Length]", false);

      if (measure) {
        improved = use_improved_estimator;

        if (!params.defined("Q0_OVER_TWO_PI")) {
          std::cerr << "parameter Q0_OVER_TWO_PI (q0/2pi) is not defined\n";
          boost::throw_exception(std::invalid_argument("correlation_length"));
        }
        if (!params.defined("DQ_OVER_TWO_PI")) {
          std::cerr << "parameter DQ_OVER_TWO_PI (dq/2pi) is not defined\n";
          boost::throw_exception(std::invalid_argument("correlation_length"));
        }

        int dim = lat.graph_helper().dimension();
        int nvs = num_sites(lat.vg());
        int nc = 2 * nvs;

        std::vector<double> q0 =
          parse_vec(static_cast<std::string>(params["Q0_OVER_TWO_PI"]), params);
        std::vector<double> dq =
          parse_vec(static_cast<std::string>(params["DQ_OVER_TWO_PI"]), params);
        if (q0.size() != dim) {
          std::cerr << "inconsistent dimension of Q0_OVER_TWO_PI (q0/2pi)\n";
          boost::throw_exception(std::invalid_argument("correlation_length::initialize"));
        }
        if (dq.size() != dim) {
          std::cerr << "inconsistent dimension of DQ_OVER_TWO_PI (dq/2pi)\n";
          boost::throw_exception(std::invalid_argument("correlation_length::initialize"));
        }
        std::clog << "info: calculating correlation length with q0 = (" << alps::write_vector(q0)
                  << ") and dq = (" << alps::write_vector(dq) << ")\n";

        phase0.resize(nvs);
        phase1.resize(nvs);
        typename alps::graph_helper<typename lattice_t::real_graph_type>::vector_type coord;
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg())) {
          coord = lat.graph_helper().coordinate(s);
          double p0 = 0;
          double p1 = 0;
          for (int i = 0; i < dim; ++i) {
            p0 += q0[i] * coord[i];
            p1 += (q0[i] + dq[i]) * coord[i];
          }
          phase0[s] = complex_type(std::cos(p0), std::sin(p0));
          phase1[s] = complex_type(std::cos(p1), std::sin(p1));
        }

        dq_abs = 0;
        for (int i = 0; i < dim; ++i) dq_abs += dq[i] * dq[i];
        dq_abs = std::sqrt(dq_abs);

        // working vector
        if (improved) {
          sq0.resize(nc);
          sq1.resize(nc);
        }

        add_scalar_obs(m, "Spin Structure Factor at q = q0", is_signed);
        add_scalar_obs(m, "Spin Structure Factor at q = q0 + dq", is_signed);
        add_scalar_obs(m, "|dq|");
      }
    }

    // improved estimator
    typedef typename dumb::template estimator<MC, LAT, TIME>::estimate estimate;
    void init_estimate(estimate&) const {}

    typedef typename dumb::template estimator<MC, LAT, TIME>::collector collector;
    void init_collector(collector&) const {}

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double /* beta */, double sign,
      std::vector<int> const& spins, std::vector<OP> const& /* operators */,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& fragments,
      collector const& /* coll */) {
      if (measure && improved) {
        int nvs = num_sites(lat.vg());
        int nc = 2 * nvs;
        std::fill(sq0.begin(), sq0.end(), complex_type(0));
        std::fill(sq1.begin(), sq1.end(), complex_type(0));
        for (int s = 0; s < nvs; ++s) {
          sq0[fragments[s].id] += phase0[s] * (0.5-spins[s]);
          sq1[fragments[s].id] += phase1[s] * (0.5-spins[s]);
        }
        double str0 = 0;
        double str1 = 0;
        for (int c = 0; c < nc; ++c) {
          str0 += power2(sq0[c]);
          str1 += power2(sq1[c]);
        }
        double vol = lat.volume();
        m["Spin Structure Factor at q = q0"] << sign * str0 / vol;
        m["Spin Structure Factor at q = q0 + dq"] << sign * str1 / vol;
        m["|dq|"] << dq_abs;
      }
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double /* beta */, double sign,
      std::vector<int> const& spins, std::vector<OP> const& /* operators */,
      std::vector<int>& /* spins_c */) {

      if (measure && !improved) {
        int nvs = num_sites(lat.vg());
        complex_type s0 = 0;
        complex_type s1 = 0;
        for (int s = 0; s < nvs; ++s) {
          s0 += phase0[s] * (0.5-spins[s]);
          s1 += phase1[s] * (0.5-spins[s]);
        }
        double vol = lat.volume();
        m["Spin Structure Factor at q = q0"] << sign * power2(s0) / vol;
        m["Spin Structure Factor at q = q0 + dq"] << sign * power2(s1) / vol;
        m["|dq|"] << dq_abs;
      }
    }
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      if (m_in.has("Spin Structure Factor at q = q0") &&
          m_in.has("Spin Structure Factor at q = q0 + dq") &&
          m_in.has("|dq|")) {
        alps::RealObsevaluator dq_abs_eval(m_in["|dq|"]);
        alps::RealObsevaluator obse_s0(m_in["Spin Structure Factor at q = q0"]);
        alps::RealObsevaluator obse_s1(m_in["Spin Structure Factor at q = q0 + dq"]);
        if (dq_abs_eval.count() && obse_s0.count() && obse_s1.count()) {
          double dq_abs = dq_abs_eval.mean();
          alps::RealObsevaluator eval("Correlation Length");
          eval = sqrt(obse_s0 / obse_s1 - 1) / dq_abs;
          m.addObservable(eval);
        }
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_CORRELATION_LENGTH_H
