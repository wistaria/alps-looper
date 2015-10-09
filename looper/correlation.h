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

#ifndef LOOPER_CORRELATION_H
#define LOOPER_CORRELATION_H

#include "measurement.h"
#include <boost/optional.hpp>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct correlation : public has_normal_estimator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<real_site_t, const typename lattice_t::virtual_graph_type,
      typename real_site_descriptor<lattice_t>::type>::type real_site_map_t;
    typedef typename alps::property_map<gauge_t, const typename lattice_t::virtual_graph_type,
      double>::type gauge_map_t;
    typedef typename alps::property_map<coordinate_t, const typename lattice_t::real_graph_type,
      coordinate_type>::type coordinate_map_t;

    bool measure_correlation, measure_green_function, measure_structure_factor, bipartite;
    real_site_map_t real_site;
    gauge_map_t gauge;
    coordinate_map_t coordinate;
    boost::optional<int> origin;
    std::valarray<double> mltplcty;
    std::vector<std::string> distance_label, momenta_label;

    // working vectors
    mutable std::valarray<double> ucorr, scorr, gucorr, gscorr, gfunc, sfac;

    void initialize(alps::Parameters const& params, lattice_t const& lat, bool /* is_signed */,
      bool enable_improved_estimator) {
      measure_correlation = params.value_or_default("MEASURE[Correlations]", false);
      measure_green_function = params.value_or_default("MEASURE[Green Function]", false);
      measure_structure_factor = params.value_or_default("MEASURE[Structure Factor]", false);
      bipartite = is_bipartite(lat);

      if (measure_green_function && !enable_improved_estimator) {
        std::cerr << "WARNING: Green funciton measurement is disabled\n";
        measure_green_function = false;
      }
      if (measure_green_function) {
        std::cerr << "WARNING: Green funciton measurement is not implemented yet\n";
        measure_green_function = false;
      }

      if (measure_correlation || measure_green_function ||
          measure_structure_factor) {
        real_site = alps::get_or_default(real_site_t(), lat.vg(), 0);
      }

      if (measure_correlation || measure_green_function) {
        gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
        if (params.defined("INITIAL_SITE"))
          origin = static_cast<int>(params["INITIAL_SITE"]);

        distance_label = distance_labels(lat, origin);
        if (!origin) {
          std::vector<unsigned int> m = distance_multiplicities(lat);
          mltplcty.resize(m.size());
          for (int i = 0; i < mltplcty.size(); ++i) mltplcty[i] = 1. / m[i];
        }
        ucorr.resize(distance_label.size());
        scorr.resize(distance_label.size());
        if (enable_improved_estimator) {
          // gucorr.resize(distance_label.size());
          // gscorr.resize(distance_label.size());
          // gfunc.resize(distance_label.size());
        }
      }

      if (measure_structure_factor) {
        coordinate = alps::get_or_default(coordinate_t(), lat.rg(), 0);
        sfac.resize(num_sites(lat.rg()));
        momenta_label = momenta_labels(lat);
      }
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      if (measure_correlation || measure_green_function) {
        if (measure_correlation) {
          add_vector_obs(m, "Spin Correlations", distance_label, is_signed);
          if (bipartite)
            add_vector_obs(m, "Staggered Spin Correlations", distance_label, is_signed);
          if (enable_improved_estimator) {
            std::cerr << "WARNING: Generalized Spin Correlations are not implemented yet\n";
            // add_vector_obs(m, "Generalized Spin Correlations", distance_label, is_signed);
            // if (bipartite)
            //   add_vector_obs(m, "Generalized Staggered Spin Correlations", distance_label,
            //     is_signed);
          }
        }
        if (measure_green_function)
          add_vector_obs(m, "Green's Function", distance_label, is_signed);
      }
      if (measure_structure_factor) {
        add_vector_obs(m, "Spin Structure Factor", momenta_label, is_signed);
      }
    }

    // improved estimator

//     typedef typename dumb::template estimator<MC, LAT, TIME>::estimate estimate;
//     void init_estimate(estimate const&) const {}

//     typedef typename dumb::template estimator<MC, LAT, TIME>::collector collector;
//     void init_collector(collector const&) const {}

//     template<typename M, typename FRAGMENT>
//     void improved_measurement(M& m,
//                               lattice_t const& lat,
//                               double /* beta */, double sign,
//                               std::vector<int> const& spins,
//                               int /* nop */,
//                               std::vector<int> const& /* spins_c */,
//                               std::vector<FRAGMENT> const& fragments,
//                               collector const& /* coll */) {
//       if ((measure_correlation || measure_green_function)) {
//         ucorr = 0;
//         scorr = 0;
//         gucorr = 0;
//         gscorr = 0;
//         gfunc = 0;
//         if (origin) {
//           typename virtual_site_iterator<lattice_t>::type si0, si0_end;
//           for (boost::tie(si0, si0_end) = sites(lat, origin.get());
//                si0 != si0_end; ++si0) {
//             double us = 1-2*spins[*si0];
//             double ss = gauge[*si0] * (1-2*spins[*si0]);
//             double gss = gauge[*si0];
//             typename virtual_site_iterator<lattice_t>::type si1, si1_end;
//             for (boost::tie(si1, si1_end) = sites(lat.vg()); si1 != si1_end;
//                  ++si1) {
//               if (fragments[*si0].id() == fragments[*si1].id()) {
//                 int r1 = real_site[*si1];
//                 ucorr[r1] += us * (1-2*spins[*si1]);
//                 scorr[r1] += ss * gauge[*si1] * (1-2*spins[*si1]);
//                 gucorr[r1] += 1;
//                 gscorr[r1] += gss * gauge[*si1];
//               }
//             }
//           }
//           ucorr *= (0.25 * sign);
//           scorr *= (0.25 * sign);
//           gucorr *= (0.25 * sign);
//           gscorr *= (0.25 * sign);
//         } else {
//           typename virtual_site_iterator<lattice_t>::type si0, si0_end;
//           for (boost::tie(si0, si0_end) = sites(lat.vg()); si0 != si0_end;
//                ++si0) {
//             int r0 = real_site[*si0];
//             double us = 1-2*spins[*si0];
//             double ss = gauge[*si0] * (1-2*spins[*si0]);
//             double gss = gauge[*si0];
//             typename virtual_site_iterator<lattice_t>::type si1, si1_end;
//             for (boost::tie(si1, si1_end) = sites(lat.vg()); si1 != si1_end;
//                  ++si1) {
//               if (fragments[*si0].id() == fragments[*si1].id()) {
//                 int d = distance(lat, r0, real_site[*si1]);
//                 ucorr[d] += us * (1-2*spins[*si1]);
//                 scorr[d] += ss * gauge[*si1] * (1-2*spins[*si1]);
//                 gucorr[d] += 1;
//                 gscorr[d] += gss * gauge[*si1];
//               }
//             }
//           }
//           ucorr *= (0.25 * sign * mltplcty);
//           scorr *= (0.25 * sign * mltplcty);
//           gucorr *= (0.25 * sign * mltplcty);
//           gscorr *= (0.25 * sign * mltplcty);
//         }
//         if (measure_correlation) {
//           m["Spin Correlations"] << ucorr;
//           m["Generalized Spin Correlations"] << gucorr;
//           if (bipartite) {
//             m["Staggered Spin Correlations"] << scorr;
//             m["Generalized Staggered Spin Correlations"] << gscorr;
//           }
//         }
//       }
//     }

//     template<typename EST, typename OP>
//     void accumulate(EST&, lattice_t const&, std::vector<int> const&,
//       std::vector<OP> const&, std::vector<int>&, double, double) {}

    struct normal_estimator {
      struct collector {
        collector() {}
        void reset(estimator_t const&) {}
        collector& operator+=(collector const&) { return *this; }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const&, lattice_t const&, double, int, int) {}
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
        template<typename M>
        void commit(M& m, estimator_t const& emt, lattice_t const& lat, double, double sign,
          double, std::vector<int> const& spins) const {
          estimator_t::real_site_map_t const& real_site = emt.real_site;
          estimator_t::gauge_map_t const& gauge = emt.gauge;
          std::valarray<double>& ucorr = emt.ucorr;
          std::valarray<double>& scorr = emt.scorr;
          std::valarray<double>& sfac = emt.sfac;
          if (emt.measure_correlation) {
            ucorr = 0;
            scorr = 0;
            if (emt.origin) {
              typename virtual_site_iterator<lattice_t>::type si0, si0_end;
              for (boost::tie(si0, si0_end) = sites(lat, emt.origin.get()); si0 != si0_end; ++si0) {
                double us = 1-2*spins[*si0];
                double ss = gauge[*si0] * (1-2*spins[*si0]);
                typename virtual_site_iterator<lattice_t>::type si1, si1_end;
                for (boost::tie(si1, si1_end) = sites(lat.vg()); si1 != si1_end; ++si1) {
                  int r1 = real_site[*si1];
                  ucorr[r1] += us * (1-2*spins[*si1]);
                  scorr[r1] += ss * gauge[*si1] * (1-2*spins[*si1]);
                }
              }
              ucorr *= 0.25 * sign;
              scorr *= 0.25 * sign;
            } else {
              typename virtual_site_iterator<lattice_t>::type si0, si0_end;
              for (boost::tie(si0, si0_end) = sites(lat.vg()); si0 != si0_end; ++si0) {
                int r0 = real_site[*si0];
                double us = 1-2*spins[*si0];
                double ss = gauge[*si0] * (1-2*spins[*si0]);
                typename virtual_site_iterator<lattice_t>::type si1, si1_end;
                for (boost::tie(si1, si1_end) = sites(lat.vg()); si1 != si1_end; ++si1) {
                  int d = distance(lat, r0, real_site[*si1]);
                  ucorr[d] += us * (1-2*spins[*si1]);
                  scorr[d] += ss * gauge[*si1] * (1-2*spins[*si1]);
                }
              }
              ucorr *= 0.25 * sign * emt.mltplcty;
              scorr *= 0.25 * sign * emt.mltplcty;
            }
            m["Spin Correlations"] << ucorr;
            if (emt.bipartite) m["Staggered Spin Correlations"] << scorr;
          }
          if (emt.measure_structure_factor) {
            sfac = 0;
            int k = 0;
            typename momentum_iterator<lattice_t>::type mit, mit_end;
            for (boost::tie(mit, mit_end) = momenta(lat); mit != mit_end; ++mit, ++k) {
              std::complex<double> val;
              BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg()))
                val += (0.5-spins[s])  * mit.phase(emt.coordinate[real_site[s]]);
              sfac[k] = sign * power2(val) / lat.volume();
            }
            m["Spin Structure Factor"] << sfac;
          }
        }
      };
    };
  };
};

} // end namespace looper

#endif // LOOPER_CORRELATION_H
