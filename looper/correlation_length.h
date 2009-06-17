/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2007-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_ONLY_PATH_INTEGRAL
# define LOOPER_ONLY_PATH_INTEGRAL
#endif

#include "measurement.h"
#include "type.h"
#include <boost/classic_spirit.hpp>

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
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef std::complex<double> complex_type;

    bool measure, improved;
    double dq_abs;
    std::vector<complex_type> phase0, phase1;

    void initialize(alps::Parameters const& params, lattice_t const& lat,
      bool /* is_signed */, bool use_improved_estimator) {
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
          coord = lat.graph_helper().coordinate(get(real_site_t(), lat.vg(), s));
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
      }
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      if (measure) {
        add_scalar_obs(m, "Spin Dynamic Structure Factor at (q,w) = (q0,0)", is_signed);
        add_scalar_obs(m, "Spin Dynamic Structure Factor at (q,w) = (q0+dq,0)", is_signed);
        add_scalar_obs(m, "|dq|");
      }
    }

    // improved estimator

    struct estimate {
      const std::vector<complex_type> *phase0_ptr;
      const std::vector<complex_type> *phase1_ptr;
      complex_type sq0, sq1;
      estimate() : sq0(complex_type(0)), sq1(complex_type(0)) {}
      void init(const std::vector<complex_type> *p0, const std::vector<complex_type> *p1) {
        phase0_ptr = p0;
        phase1_ptr = p1;
        sq0 = complex_type(0);
        sq1 = complex_type(0);
      }
      estimate& operator+=(estimate const& rhs) {
        sq0 += rhs.sq0;
        sq1 += rhs.sq1;
        return *this;
      }
      void begin_s(lattice_t const& lat, double t, int s, int c) { end_s(lat, -t, s, c); }
      void begin_bs(lattice_t const& lat, double t, int, int s, int c) { begin_s(lat, t, s, c); }
      void begin_bt(lattice_t const& lat, double t, int, int s, int c) { begin_s(lat, t, s, c); }
      void end_s(lattice_t const&, double t, int s, int c) {
        sq0 += (*phase0_ptr)[s] * t * (0.5-c);
        sq1 += (*phase1_ptr)[s] * t * (0.5-c);
      }
      void end_bs(lattice_t const& lat, double t, int, int s, int c) { end_s(lat, t, s, c); }
      void end_bt(lattice_t const& lat, double t, int, int s, int c) { end_s(lat, t, s, c); }
      void start_bottom(lattice_t const& lat, double t, int s, int c) { begin_s(lat, t, s, c); }
      void start(lattice_t const& lat, double t, int s, int c) { begin_s(lat, t, s, c); }
      void stop(lattice_t const& lat, double t, int s, int c) { end_s(lat, t, s, c); }
      void stop_top(lattice_t const& lat, double t, int s, int c) { end_s(lat, t, s, c); }
    };
    void init_estimate(estimate& est) const { est.init(&phase0, &phase1); }

    struct collector {
      double str0, str1;
      void init() {
        str0 = 0; str1 = 0;
      }
      collector& operator+=(collector const& coll) {
        str0 += coll.str0;
        str1 += coll.str1;
        return *this;
      }
      collector& operator+=(estimate const& est) {
        str0 += power2(est.sq0);
        str1 += power2(est.sq1);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const& lat, double beta, int /* nop */, double sign) const {
        double vol = lat.volume();
        m["Spin Dynamic Structure Factor at (q,w) = (q0,0)"] << sign * beta * str0 / vol;
        m["Spin Dynamic Structure Factor at (q,w) = (q0+dq,0)"] << sign * beta * str1 / vol;
      }
    };
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& /* spins */, std::vector<OP> const& operators,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
      m["|dq|"] << dq_abs;
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators,
      std::vector<int>& spins_c) {
      if (!measure || improved) return;

      complex_type sq0 = 0;
      complex_type sq1 = 0;
      typename virtual_site_iterator<lattice_t>::type si, si_end;
      for (boost::tie(si, si_end) = sites(lat.vg()); si != si_end; ++si) {
        sq0 += phase0[*si] * (0.5-spins[*si]);
        sq1 += phase1[*si] * (0.5-spins[*si]);
      }
      complex_type sq0_a = 0;
      complex_type sq1_a = 0;
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          double t = oi->time();
          sq0_a += t * sq0;
          sq1_a += t * sq1;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            sq0 += phase0[s] * (1.0-2*spins_c[s]);
            sq1 += phase1[s] * (1.0-2*spins_c[s]);
          } else {
            unsigned int s0 = source(oi->pos(), lat.vg());
            unsigned int s1 = target(oi->pos(), lat.vg());
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            sq0 += phase0[s0] * (1.0-2*spins_c[s0]) + phase0[s1] * (1.0-2*spins_c[s1]);
            sq1 += phase1[s0] * (1.0-2*spins_c[s0]) + phase1[s1] * (1.0-2*spins_c[s1]);
          }
          sq0_a -= t * sq0;
          sq1_a -= t * sq1;
        }
      }
      sq0_a += sq0;
      sq1_a += sq1;
      double vol = lat.volume();
      m["Spin Dynamic Structure Factor at (q,w) = (q0,0)"]
        << sign * beta * power2(sq0_a) / vol;
      m["Spin Dynamic Structure Factor at (q,w) = (q0+dq,0)"]
        << sign * beta * power2(sq1_a) / vol;
      m["|dq|"] << dq_abs;
    }
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      if (m_in.has("Spin Dynamic Structure Factor at (q,w) = (q0,0)") &&
          m_in.has("Spin Dynamic Structure Factor at (q,w) = (q0+dq,0)") &&
          m_in.has("|dq|")) {
        alps::RealObsevaluator dq_abs_eval(m_in["|dq|"]);
        alps::RealObsevaluator obse_s0(m_in["Spin Dynamic Structure Factor at (q,w) = (q0,0)"]);
        alps::RealObsevaluator obse_s1(m_in["Spin Dynamic Structure Factor at (q,w) = (q0+dq,0)"]);
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
