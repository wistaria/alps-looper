/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TWIST_ORDER_PARAMETER_H
#define LOOPER_TWIST_ORDER_PARAMETER_H

#include "measurement.h"
#include <boost/foreach.hpp>
#include <complex>

namespace looper {

struct twist_order_parameter {

  typedef dumb_measurement<twist_order_parameter> dumb;
  
  template<typename MC, typename LAT, typename TIME>
  struct estimator
  {
    typedef MC   mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<gauge_t,
              const typename lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool improved;
    std::vector<double> phase;

    template<typename M>
    void initialize(M& m, alps::Parameters const& /* params */,
                    lattice_t const& lat,
                    bool is_signed, bool use_improved_estimator)
    {
      improved = use_improved_estimator;

      // check basis vector
      int i = 0;
      BOOST_FOREACH(std::vector<double> const& v, lat.graph_helper().basis_vectors()) {
        int j = 0;
        BOOST_FOREACH(double x, v) {
          if ((j == i && alps::is_zero(x)) || (j != i && alps::is_nonzero(x)))
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
      
      add_scalar_obs(m, "Twist Order Parameter", is_signed);
      add_scalar_obs(m, "Twist Order Parameter (Imaginary Part)", is_signed);
    }

    // improved estimator

    struct estimate {
      const std::vector<double> *phase_ptr;
      double moment;
      void init(const std::vector<double> *p) {
        phase_ptr = p;
        moment = 0;
      }
      void start_s(lattice_t const&, double, int, int) {}
      void start_bs(lattice_t const&, double, int, int, int) {}
      void start_bt(lattice_t const&, double, int, int, int) {}
      void term_s(lattice_t const&, double, int, int) {}
      void term_bs(lattice_t const&, double, int, int, int) {}
      void term_bt(lattice_t const&, double, int, int, int) {}
      void at_bot(lattice_t const&, double, int s, int c) {
        moment += (*phase_ptr)[s] * (0.5-c);
      }
      void at_top(lattice_t const&, double, int, int) {}
    };
    void init_estimate(estimate& est) const { est.init(&phase); }
      
    struct collector {
      double top;
      void init() { top = 1; }
      template<typename EST>
      collector operator+(EST const& est) {
        top *= std::cos(est.moment);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const&, double, int, double sign) const {
        m["Twist Order Parameter"] << sign * top;
        m["Twist Order Parameter (Imaginary Part)"] << 0.;
      }
    };
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m,
                              lattice_t const& lat,
                              double beta, double sign,
                              std::vector<int> const& /* spins */,
                              std::vector<OP> const& operators,
                              std::vector<int> const& /* spins_c */,
                              std::vector<FRAGMENT> const& /* fragments */,
                              collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat,
                            double /* beta */, double sign,
                            std::vector<int> const& spins,
                            std::vector<OP> const& /* operators */,
                            std::vector<int>& /* spins_c */)
    {
      if (improved) return;

      double total = 0;
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg()))
        total += phase[s] * (0.5-spins[s]);
      m["Twist Order Parameter"] << sign * std::cos(total);
      m["Twist Order Parameter (Imaginary Part)"] << sign * std::sin(total);
    }
  };

  typedef dumb::evaluator evaluator;
};

} // end namespace looper

#endif // LOOPER_TWIST_ORDER_PARAMTER_H
