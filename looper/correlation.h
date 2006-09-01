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

#ifndef LOOPER_CORRELATION_H
#define LOOPER_CORRELATION_H

#include "measurement.h"
#include <boost/optional.hpp>

namespace looper {

struct correlation
{
  typedef dumb_measurement<correlation> dumb;

  template<typename MC, typename VLAT, typename TIME>
  struct estimator
  {
    typedef MC   mc_type;
    typedef VLAT virtual_lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<gauge_t,
              const typename virtual_lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool measure;
    bool improved;
    gauge_map_t gauge;
    boost::optional<typename virtual_lattice_t::real_site_descriptor> origin;
    std::valarray<double> mltplcty;

    // working vector
    std::valarray<double> ucorr, scorr, gucorr, gscorr;

    template<typename M>
    void initialize(M& m, alps::Parameters const& params,
                    virtual_lattice_t const& vlat,
                    bool is_signed, bool use_improved_estimator)
    {
      measure = params.value_or_default("MEASURE[Correlations]", false);
      if (!measure) return;

      improved = use_improved_estimator;
      gauge = alps::get_or_default(gauge_t(), vlat.vgraph(), 0);
      if (params.defined("INITIAL_SITE"))
        origin = static_cast<int>(params["INITIAL_SITE"]);

      std::vector<std::string> label;
      if (origin) {
        for (int s = 0; s < num_vertices(vlat.rgraph()); ++s)
          label.push_back(vlat.helper().coordinate_string(origin.get()) +
                          " -- " + vlat.helper().coordinate_string(s));
        ucorr.resize(vlat.helper().num_sites());
        scorr.resize(vlat.helper().num_sites());
        if (improved) {
          gucorr.resize(vlat.helper().num_sites());
          gscorr.resize(vlat.helper().num_sites());
        }
      } else {
        label = vlat.helper().distance_labels();
        std::vector<unsigned int> m = vlat.helper().distance_multiplicities();
        mltplcty.resize(m.size());
        for (int i = 0; i < mltplcty.size(); ++i) mltplcty[i] = 1.0 / m[i];
        ucorr.resize(vlat.helper().num_distances());
        scorr.resize(vlat.helper().num_distances());
        if (improved) {
          gucorr.resize(vlat.helper().num_distances());
          gscorr.resize(vlat.helper().num_distances());
        }
      }

      add_vector_obs(m, "Spin Correlations", label, is_signed);
      if (improved)
        add_vector_obs(m, "Generalized Spin Correlations", label, is_signed);
      if (is_bipartite(vlat)) {
        add_vector_obs(m, "Staggered Spin Correlations", label, is_signed);
        if (improved)
          add_vector_obs(m, "Generalized Staggered Spin Correlations",
                         label, is_signed);
      }
    }

    // improved estimator

    typedef typename dumb::template estimator<MC, VLAT, TIME>::estimate
      estimate;
    void init_estimate(estimate const&) const {}

    typedef typename dumb::template estimator<MC, VLAT, TIME>::collector
      collector;
    void init_collector(collector const&) const {}

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m,
                              virtual_lattice_t const& vlat,
                              double /* beta */, double sign,
                              std::vector<int> const& spins,
                              std::vector<OP> const& /* operators */,
                              std::vector<int> const& /* spins_c */,
                              std::vector<FRAGMENT> const& fragments,
                              collector const& /* coll */)
    {
      if (!measure || !improved) return;

      ucorr = 0;
      scorr = 0;
      gucorr = 0;
      gscorr = 0;
      if (origin) {
        typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
        for (boost::tie(si0, si0_end) = virtual_sites(vlat, origin.get());
             si0 != si0_end; ++si0) {
          double us = 0.25 * (1-2*spins[*si0]);
          double ss = 0.25 * gauge[*si0] * (1-2*spins[*si0]);
          double gss = 0.25 * gauge[*si0];
          typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
          for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end; ++si1) {
            if (fragments[*si0].id == fragments[*si1].id) {
              int r1 = rsite(vlat, *si1);
              ucorr[r1] += us * (1-2*spins[*si1]);
              scorr[r1] += ss * gauge[*si1] * (1-2*spins[*si1]);
              gucorr[r1] += 0.25;
              gscorr[r1] += gss * gauge[*si1];
            }
          }
        }
        m["Spin Correlations"] << std::valarray<double>(sign * ucorr);
        m["Staggered Spin Correlations"]
          << std::valarray<double>(sign * scorr);
        m["Generalized Spin Correlations"]
          << std::valarray<double>(sign * gucorr);
        m["Generalized Staggered Spin Correlations"]
          << std::valarray<double>(sign * gscorr);
      } else {
        typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
        for (boost::tie(si0, si0_end) = vsites(vlat); si0 != si0_end; ++si0) {
          int r0 = rsite(vlat, *si0);
          double us = 0.25 * (1-2*spins[*si0]);
          double ss = 0.25 * gauge[*si0] * (1-2*spins[*si0]);
          double gss = 0.25 * gauge[*si0];
          typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
          for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end; ++si1) {
            if (fragments[*si0].id == fragments[*si1].id) {
              int d = vlat.helper().distance(r0, rsite(vlat, *si1));
              ucorr[d] += us * (1-2*spins[*si1]);
              scorr[d] += ss * gauge[*si1] * (1-2*spins[*si1]);
              gucorr[d] += 0.25;
              gscorr[d] += gss * gauge[*si1];
            }
          }
        }
        m["Spin Correlations"]
          << std::valarray<double>(sign * mltplcty * ucorr);
        m["Staggered Spin Correlations"]
          << std::valarray<double>(sign * mltplcty * scorr);
        m["Generalized Spin Correlations"]
          << std::valarray<double>(sign * mltplcty * gucorr);
        m["Generalized Staggered Spin Correlations"]
          << std::valarray<double>(sign * mltplcty * gscorr);
      }
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, virtual_lattice_t const& vlat,
                            double /* beta */, double sign,
                            std::vector<int> const& spins,
                            std::vector<OP> const& /* operators */,
                            std::vector<int>& /* spins_c */)
    {
      if (!measure || improved) return;

      ucorr = 0;
      scorr = 0;
      if (origin) {
        typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
        for (boost::tie(si0, si0_end) = virtual_sites(vlat, origin.get());
             si0 != si0_end; ++si0) {
          double us = 0.25 * (1-2*spins[*si0]);
          double ss = 0.25 * gauge[*si0] * (1-2*spins[*si0]);
          typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
          for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end; ++si1) {
            int r1 = rsite(vlat, *si1);
            ucorr[r1] += us * (1-2*spins[*si1]);
            scorr[r1] += ss * gauge[*si1] * (1-2*spins[*si1]);
          }
        }
        m["Spin Correlations"] << std::valarray<double>(sign * ucorr);
        m["Staggered Correlations"] << std::valarray<double>(sign * scorr);
      } else {
        typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
        for (boost::tie(si0, si0_end) = vsites(vlat); si0 != si0_end; ++si0) {
          int r0 = rsite(vlat, *si0);
          double us = 0.25 * (1-2*spins[*si0]);
          double ss = 0.25 * gauge[*si0] * (1-2*spins[*si0]);
          typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
          for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end; ++si1) {
            int d = vlat.helper().distance(r0, rsite(vlat, *si1));
            ucorr[d] += us * (1-2*spins[*si1]);
            scorr[d] += ss * gauge[*si1] * (1-2*spins[*si1]);
          }
        }
        m["Spin Correlations"]
          << std::valarray<double>(sign * mltplcty * ucorr);
        m["Staggered Correlations"]
          << std::valarray<double>(sign * mltplcty * scorr);
      }
    }
  };

  typedef dumb::evaluator evaluator;
};

} // end namespace looper

#endif // LOOPER_CORRELATION_H
