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
    typedef typename alps::property_map<alps::coordinate_t,
              const typename virtual_lattice_t::real_graph_type,
              alps::coordinate_type>::type coordinate_map_t;
    typedef typename alps::property_map<gauge_t,
              const typename virtual_lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool measure_correlation, measure_green_function,
      measure_structure_factor;
    bool improved;
    gauge_map_t gauge;
    coordinate_map_t coordinate;
    boost::optional<typename virtual_lattice_t::real_site_descriptor> origin;
    std::valarray<double> mltplcty;

    // working vector
    std::valarray<double> ucorr, scorr, gucorr, gscorr, gfunc, sfac;

    template<typename M>
    void initialize(M& m, alps::Parameters const& params,
                    virtual_lattice_t const& vlat,
                    bool is_signed, bool use_improved_estimator)
    {
      measure_correlation =
        params.value_or_default("MEASURE[Correlations]", false);
      measure_green_function =
        params.value_or_default("MEASURE[Green Function]", false);
      measure_structure_factor =
        params.value_or_default("MEASURE[Structure Factor]", false);
      improved = use_improved_estimator;

      if (measure_green_function && !improved) {
        std::cerr << "Warning: Green funciton measurement is disabled\n";
        measure_green_function = false;
      }

      if (measure_correlation || measure_green_function) {
        gauge = alps::get_or_default(gauge_t(), vlat.vgraph(), 0);
        if (params.defined("INITIAL_SITE"))
          origin = static_cast<int>(params["INITIAL_SITE"]);

        std::vector<std::string> label;
        if (origin) {
          for (int s = 0; s < num_vertices(vlat.rgraph()); ++s)
            label.push_back(vlat.helper().coordinate_string(origin.get()) +
                            " -- " + vlat.helper().coordinate_string(s));
        } else {
          label = vlat.helper().distance_labels();
          std::vector<unsigned int> m =
            vlat.helper().distance_multiplicities();
          mltplcty.resize(m.size());
          for (int i = 0; i < mltplcty.size(); ++i) mltplcty[i] = 1.0 / m[i];
        }
        if (measure_correlation || measure_green_function) {
          ucorr.resize(label.size());
          scorr.resize(label.size());
          if (improved) {
            gucorr.resize(label.size());
            gscorr.resize(label.size());
            gfunc.resize(label.size());
          }
        }
        if (measure_correlation) {
          add_vector_obs(m, "Spin Correlations", label, is_signed);
          if (is_bipartite(vlat))
            add_vector_obs(m, "Staggered Spin Correlations", label, is_signed);
          if (improved) {
            add_vector_obs(m, "Generalized Spin Correlations", label,
                           is_signed);
            if (is_bipartite(vlat))
              add_vector_obs(m, "Generalized Staggered Spin Correlations",
                             label, is_signed);
          }
        }
        if (measure_green_function)
          add_vector_obs(m, "Green's Function", label, is_signed);
      }

      if (measure_structure_factor) {
        coordinate =
          alps::get_or_default(alps::coordinate_t(), vlat.helper().graph(), 0);
        sfac.resize(vlat.helper().num_sites());
        add_vector_obs(m, "Spin Structure Factor",
                       vlat.helper().momenta_labels(), is_signed);
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
      if ((measure_correlation || measure_green_function) && improved) {
        ucorr = 0;
        scorr = 0;
        gucorr = 0;
        gscorr = 0;
        gfunc = 0;
        if (origin) {
          typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
          for (boost::tie(si0, si0_end) = virtual_sites(vlat, origin.get());
               si0 != si0_end; ++si0) {
            double us = 1-2*spins[*si0];
            double ss = gauge[*si0] * (1-2*spins[*si0]);
            double gss = gauge[*si0];
            typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
            for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end;
                 ++si1) {
              if (fragments[*si0].id == fragments[*si1].id) {
                int r1 = rsite(vlat, *si1);
                ucorr[r1] += us * (1-2*spins[*si1]);
                scorr[r1] += ss * gauge[*si1] * (1-2*spins[*si1]);
                gucorr[r1] += 1;
                gscorr[r1] += gss * gauge[*si1];
              }
            }
          }
          ucorr *= (0.25 * sign);
          scorr *= (0.25 * sign);
          gucorr *= (0.25 * sign);
          gscorr *= (0.25 * sign);
        } else {
          typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
          for (boost::tie(si0, si0_end) = vsites(vlat); si0 != si0_end;
               ++si0) {
            int r0 = rsite(vlat, *si0);
            double us = 1-2*spins[*si0];
            double ss = gauge[*si0] * (1-2*spins[*si0]);
            double gss = gauge[*si0];
            typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
            for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end;
                 ++si1) {
              if (fragments[*si0].id == fragments[*si1].id) {
                int d = vlat.helper().distance(r0, rsite(vlat, *si1));
                ucorr[d] += us * (1-2*spins[*si1]);
                scorr[d] += ss * gauge[*si1] * (1-2*spins[*si1]);
                gucorr[d] += 1;
                gscorr[d] += gss * gauge[*si1];
              }
            }
          }
          ucorr *= (0.25 * sign * mltplcty);
          scorr *= (0.25 * sign * mltplcty);
          gucorr *= (0.25 * sign * mltplcty);
          gscorr *= (0.25 * sign * mltplcty);
        }
        if (measure_correlation) {
          m["Spin Correlations"] << ucorr;
          m["Generalized Spin Correlations"] << gucorr;
          if (is_bipartite(vlat)) {
            m["Staggered Spin Correlations"] << scorr;
            m["Generalized Staggered Spin Correlations"] << gscorr;
          }
        }
      }
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, virtual_lattice_t const& vlat,
                            double /* beta */, double sign,
                            std::vector<int> const& spins,
                            std::vector<OP> const& /* operators */,
                            std::vector<int>& /* spins_c */)
    {
      if (measure_correlation && !improved) {
        ucorr = 0;
        scorr = 0;
        if (origin) {
          typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
          for (boost::tie(si0, si0_end) = virtual_sites(vlat, origin.get());
               si0 != si0_end; ++si0) {
            double us = 1-2*spins[*si0];
            double ss = gauge[*si0] * (1-2*spins[*si0]);
            typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
            for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end;
                 ++si1) {
              int r1 = rsite(vlat, *si1);
              ucorr[r1] += us * (1-2*spins[*si1]);
              scorr[r1] += ss * gauge[*si1] * (1-2*spins[*si1]);
            }
          }
          ucorr *= 0.25 * sign;
          scorr *= 0.25 * sign;
        } else {
          typename virtual_lattice_t::virtual_site_iterator si0, si0_end;
          for (boost::tie(si0, si0_end) = vsites(vlat); si0 != si0_end;
               ++si0) {
            int r0 = rsite(vlat, *si0);
            double us = 1-2*spins[*si0];
            double ss = gauge[*si0] * (1-2*spins[*si0]);
            typename virtual_lattice_t::virtual_site_iterator si1, si1_end;
            for (boost::tie(si1, si1_end) = vsites(vlat); si1 != si1_end;
                 ++si1) {
              int d = vlat.helper().distance(r0, rsite(vlat, *si1));
              ucorr[d] += us * (1-2*spins[*si1]);
              scorr[d] += ss * gauge[*si1] * (1-2*spins[*si1]);
            }
          }
          ucorr *= 0.25 * sign * mltplcty;
          scorr *= 0.25 * sign * mltplcty;
        }
        m["Spin Correlations"] << ucorr;
        if (is_bipartite(vlat))
          m["Staggered Spin Correlations"] << scorr;
      }

      if (measure_structure_factor) {
        sfac = 0;
        int k = 0;
        typename virtual_lattice_t::graph_helper_type::momentum_iterator
          mit, mit_end;
        for (boost::tie(mit, mit_end) = vlat.helper().momenta();
             mit != mit_end; ++mit, ++k) {
          std::complex<double> val;
          typename virtual_lattice_t::virtual_site_iterator si, si_end;
          for (boost::tie(si, si_end) = vvertices(vlat); si != si_end; ++si)
            val += (0.5-spins[*si])  * mit.phase(coordinate[rsite(vlat, *si)]);
          sfac[k] = sign * power2(val) / vlat.helper().num_sites();
        }
        m["Spin Structure Factor"] << sfac;
      }
    }
  };

  typedef dumb::evaluator evaluator;
};

} // end namespace looper

#endif // LOOPER_CORRELATION_H
