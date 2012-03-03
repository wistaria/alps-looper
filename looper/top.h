/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TOP_H
#define LOOPER_TOP_H

#include "measurement.h"
#include <alps/numeric/is_nonzero.hpp>
#include <alps/numeric/is_zero.hpp>
#include <boost/foreach.hpp>
#include <complex>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

template<unsigned int N>
struct twist_order_parameter_n : public has_improved_estimator_tag,
  public has_normal_estimator_tag {

  BOOST_STATIC_ASSERT(N > 0);

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;

    std::vector<double> phase;
    std::vector<std::string> label;

    void initialize(alps::Parameters const& /* params */, lattice_t const& lat,
      bool /* is_signed */, bool /* enable_improved_estimator */) {
      // check basis vector
      int i = 0;
      BOOST_FOREACH(std::vector<double> const& v, lat.graph_helper().basis_vectors()) {
        int j = 0;
        BOOST_FOREACH(double x, v) {
          if ((j == i && is_zero(x)) || (j != i && is_nonzero(x)))
            boost::throw_exception(std::runtime_error("basis vector check failed"));
          ++j;
        }
        ++i;
      }

      const double pi = std::acos(-1.);
      const double span =
        lat.graph_helper().lattice().extent(0) * lat.graph_helper().basis_vectors().first->front();
      phase.clear();
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg()))
        phase.push_back(2 * pi * get(coordinate_t(), lat.rg(),
          get(real_site_t(), lat.vg(), s)).front() / span);

      // std::cerr << "span is: " << span << std::endl;
      // BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg()))
      //   std::cerr << s << ' ' << phase[s] << std::endl;

      label.resize(N);
      label[0] = "Twist Order Parameter";
      for (unsigned int i = 1; i < N; ++i)
        label[i] = "Twist Order Parameter with p = " + boost::lexical_cast<std::string>(i+1);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool /* enable_improved_estimator */) {
      for (unsigned int i = 0; i < label.size(); ++i) {
        add_scalar_obs(m, label[i], is_signed);
        add_scalar_obs(m, label[i] + " (Imaginary Part)", is_signed);
      }
    }

    struct improved_estimator {
      struct estimate {
        double moment;
        estimate() : moment(0) {}
        void reset(estimator_t const&) { moment = 0; }
        estimate& operator+=(estimate const& rhs) { moment += rhs.moment; }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const& emt, lattice_t const&, double, int s, int c) {
          moment += emt.phase[s] * (0.5-c);
        }
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
      };

      struct collector {
        double top;
        void reset(estimator_t const&) { top = 1; }
        collector& operator+=(collector const& coll) {
          top *= coll.top;
          return *this;
        }
        collector& operator+=(estimate const& est) {
          top *= std::cos(est.moment);
          return *this;
        }
        template<typename M>
        void commit(M& m, estimator_t const& emt, lattice_t const&, double, double sign,
          double, std::vector<int> const&) const {
          double t = sign * top;
          for (int i = 0; i < N; ++i) {
            m[emt.label[i]] << t;
            m[emt.label[i] + " (Imaginary Part)"] << 0.;
            t *= top;
          }
        }
      };
    };

    struct normal_estimator {
      struct collector {
        double moment;
        collector() : moment(0) {}
        void reset(estimator_t const&) { moment = 0; }
        collector& operator+=(collector const& rhs) {
          moment += rhs.moment;
          return *this;
        }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const& emt, lattice_t const&, double, int s, int c) {
          moment += emt.phase[s] * (0.5-c);
        }
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
        template<typename M>
        void commit(M& m, estimator_t const& emt, lattice_t const&, double, double sign,
          double, std::vector<int> const&) const {
          for (int i = 0; i < N; ++i) {
            m[emt.label[i]] << sign * std::cos((i+1) * moment);
            m[emt.label[i] + " (Imaginary Part)"] << sign * std::sin((i+1) * moment);
          }
        }
      };
    };
  };
};

typedef twist_order_parameter_n<1> twist_order_parameter;

} // end namespace looper

#endif // LOOPER_TOP_H
