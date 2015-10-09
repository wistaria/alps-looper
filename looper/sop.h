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

#ifndef LOOPER_SOP_H
#define LOOPER_SOP_H

#include "measurement.h"
#include <alps/numeric/is_nonzero.hpp>
#include <alps/numeric/is_zero.hpp>
#include <boost/foreach.hpp>
#include <cmath>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct string_order_parameter : public has_normal_estimator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;

    int p_left;
    int p_right;
    std::string label;

    void initialize(alps::Parameters const& params, lattice_t const& lat,
      bool /* is_signed */, bool /* enable_improved_estimator */) {
      p_left = 1;
      p_right = 1;
      if (params.defined("SOP_M") && params.defined("SOP_N")) {
        p_left = static_cast<int>(evaluate(params["SOP_M"], params));
        p_right = static_cast<int>(evaluate(params["SOP_N"], params));
      }

      // check basis vector
      int i = 0;
      BOOST_FOREACH(std::vector<double> const& v, lat.graph_helper().basis_vectors()) {
        int j = 0;
        BOOST_FOREACH(double x, v) {
          if ((j == i && alps::numeric::is_zero(x)) || (j != i && alps::numeric::is_nonzero(x)))
            boost::throw_exception(std::runtime_error("basis vector check failed"));
          ++j;
        }
        ++i;
      }

      // check parity
      if (num_sites(lat.rg()) < 4)
        boost::throw_exception(std::runtime_error("sop: too small lattice"));
      if (num_sites(lat.rg()) % 1 != 0)
        boost::throw_exception(std::runtime_error("sop: number of sites should be odd"));
      int s2 = 0;
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, 0)) ++s2;
      for (int r = 1; r < num_sites(lat.rg()); ++r) {
        int n = 0;
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, r)) ++n;
        if (n != s2)
          boost::throw_exception(std::runtime_error("sop: only uniform spins supported"));
      }
      if ((s2 + p_left + p_right) & 1 != 0)
          boost::throw_exception(std::runtime_error("sop: parity error"));

      label = "String Order Parameter";
      if (params.defined("SOP_M") && params.defined("SOP_N"))
        label += " (M,N)=(" + boost::lexical_cast<std::string>(p_left) + "," +
          boost::lexical_cast<std::string>(p_right) + ")";
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool /* enable_improved_estimator */) {
      add_scalar_obs(m, label, is_signed);
    }

    struct normal_estimator {
      struct collector {
        collector() {}
        void reset(estimator_t const&) {}
        collector& operator+=(collector const& rhs) { return *this; }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const& emt, lattice_t const&, double, int s, int c) {}
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
        template<typename M>
        void commit(M& m, estimator_t const& emt, lattice_t const& lat, double /* beta */, double sign, double,
          std::vector<int> const& spins) const {
          using std::cos; using std::pow;
          int left = 0;
          int right = num_sites(lat.rg()) / 2;
          int regsz2 = 0; // sum of 2Sz in the region between left and right
          for (int r = left; r < right; ++r)
            BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, r))
              regsz2 += 1 - 2 * spins[v];

          double sum = 0;
          for (int s = 0; s < num_sites(lat.rg()); ++s) {
            int lsz2 = 0;
            BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, left))
              lsz2 += 1 - 2 * spins[v];
            int rsz2 = 0;
            BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, right))
              rsz2 += 1 - 2 * spins[v];
            regsz2 -= lsz2;
            if ((left & 1) == 0) {
              sum += pow(0.5 * lsz2, (double)emt.p_left) *
                cos(0.5 * M_PI * (regsz2 + emt.p_left + emt.p_right)) *
                pow(0.5 * rsz2, (double)emt.p_right);
            }
            regsz2 += rsz2;
            ++left;
            ++right;
            if (right == num_sites(lat.rg())) right = 0;
          }
          m[emt.label] << sign * sum / (num_sites(lat.rg()) / 2);
        }
      };
    };
  };
};

} // end namespace looper

#endif // LOOPER_SOP_H
