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

#ifndef LOOPER_SOP_H
#define LOOPER_SOP_H

#include "measurement.h"
#include <boost/foreach.hpp>
#include <cmath>

namespace looper {

struct string_order_parameter {

  typedef dumb_measurement<string_order_parameter> dumb;

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;

    int p_left;
    int p_right;
    std::string label;

    template<typename M>
    void initialize(M& m, alps::Parameters const& params, lattice_t const& lat,
      bool is_signed, bool /* use_improved_estimator */) {

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

      // check parity
      if (num_sites(lat.rg()) < 4)
        boost::throw_exception(std::runtime_error("sop: too small lattice"));
      int left = 0;
      int right = num_sites(lat.rg()) / 2 - 1;
      int ns = 0; // number of spins between left and right
      for (int r = left; r < right; ++r)
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, r)) ++ns;
      for (int s = 0; s < num_sites(lat.rg()); ++s) {
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, left)) --ns;
        if (ns & 1 == 1) boost::throw_exception(std::runtime_error("sop: parity error"));
        ++left;
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, right)) ++ns;
        ++right;
        if (right == num_sites(lat.rg())) right = 0;
      }

      p_left = 1;
      p_right = 1;
      label = "String Order Parameter";
      if (params.defined("SOP_M") && params.defined("SOP_N")) {
        p_left = static_cast<int>(evaluate(params["SOP_M"], params));
        p_right = static_cast<int>(evaluate(params["SOP_N"], params));
        label += " (M,N)=(" + boost::lexical_cast<std::string>(p_left) + "," +
          boost::lexical_cast<std::string>(p_right) + ")";
      }
      add_scalar_obs(m, label, is_signed);
    }

    // improved estimator

    typedef typename dumb::template estimator<MC, LAT, TIME>::estimate estimate;
    void init_estimate(estimate&) const {}

    typedef typename dumb::template estimator<MC, LAT, TIME>::collector collector;
    void init_collector(collector&) const {}

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M&, lattice_t const&, double, double, std::vector<int> const&,
      std::vector<OP> const&, std::vector<int> const&, std::vector<FRAGMENT> const&,
      collector const&) {}

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double /* beta */, double sign,
      std::vector<int> const& spins, std::vector<OP> const& /* operators */,
      std::vector<int>& /* spins_c */) {

      using std::cos; using std::pow;

      int left = 0;
      int right = num_sites(lat.rg()) / 2 - 1;
      int regsz2 = 0; // sum of 2Sz in the region between left and right
      for (int r = left; r < right; ++r)
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, r))
          regsz2 += 1 - 2 * spins[v];

      double sum = 0;
      for (int s = 0; s < num_sites(lat.rg()); ++s) {
        int lsz2 = 0;
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, left)) {
          regsz2 -= 1 - 2 * spins[v];
          lsz2 += 1 - 2 * spins[v];
        }
        int rsz2 = 0;
        BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type v, sites(lat, right))
          rsz2 += 1 - 2 * spins[v];
        sum += pow(0.5 * lsz2, (double)p_left) * cos(0.5 * M_PI * regsz2) * pow(0.5 * rsz2, (double)p_right);
        ++left;
        regsz2 += rsz2;
        ++right;
        if (right == num_sites(lat.rg())) right = 0;
      }
      m[label] << sign * sum / num_sites(lat.rg());
    }
  };

  typedef dumb::evaluator evaluator;
};

} // end namespace looper

#endif // LOOPER_SOP_H
