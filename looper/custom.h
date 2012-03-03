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

#ifndef LOOPER_CUSTOM_H
#define LOOPER_CUSTOM_H

#include "measurement.h"
#include <alps/scheduler/measurement_operators.h>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <vector>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

template<typename LAT>
struct custom_measurement_initializer : private alps::MeasurementOperators {
  typedef LAT lattice_t;
  typedef boost::tuple<std::string, std::vector<std::vector<double> > > s_elements_type;
  typedef boost::tuple<std::string, std::vector<std::vector<double> >,
    std::vector<std::vector<double> > > p_elements_type;
  typedef std::pair<std::string, std::string> s_expression_type;
  typedef std::pair<std::string, std::pair<std::string, std::string> > p_expression_type;

  custom_measurement_initializer(alps::Parameters const& params);

  void init(lattice_t const& lat,
    std::vector<s_elements_type>& average_elements,
    std::vector<s_elements_type>& local_elements,
    std::vector<p_elements_type>& correlation_elements,
    std::vector<p_elements_type>& strfactor_elements);

private:
  alps::Parameters const& params_;
};

struct custom_measurement : has_normal_estimator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename custom_measurement_initializer<lattice_t>::s_elements_type s_elements_type;
    typedef typename custom_measurement_initializer<lattice_t>::p_elements_type p_elements_type;
    typedef typename alps::property_map<coordinate_t, const typename lattice_t::real_graph_type,
      coordinate_type>::type coordinate_map_t;

    std::vector<s_elements_type> average_elements;
    std::vector<s_elements_type> local_elements;
    std::vector<p_elements_type> correlation_elements;
    std::vector<p_elements_type> strfactor_elements;
    coordinate_map_t coordinate;
    boost::optional<int> origin;
    std::valarray<double> mltplcty;
    std::vector<std::string> site_label, distance_label, momenta_label;

    // working vectors
    mutable std::valarray<double> local, corr, sfac;

    void initialize(alps::Parameters const& params, lattice_t const& lat, bool /* is_signed */,
      bool /* enable_improved_estimator */) {
      custom_measurement_initializer<lattice_t> initializer(params);
      initializer.init(lat, average_elements, local_elements, correlation_elements,
                       strfactor_elements);
      if (local_elements.size()) {
        site_label = alps::site_labels(lat.rg());
        local.resize(site_label.size());
      }
      if (correlation_elements.size()) {
        if (params.defined("INITIAL_SITE"))
          origin = static_cast<int>(params["INITIAL_SITE"]);
        distance_label = distance_labels(lat, origin);
        if (!origin) {
          std::vector<unsigned int> m = distance_multiplicities(lat);
          mltplcty.resize(m.size());
          for (int i = 0; i < mltplcty.size(); ++i) mltplcty[i] = 1. / m[i];
        }
        corr.resize(distance_label.size());
      }
      if (strfactor_elements.size()) {
        coordinate = alps::get_or_default(coordinate_t(), lat.rg(), 0);
        momenta_label = momenta_labels(lat);
        sfac.resize(momenta_label.size());
      }
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool /* enable_improved_estimator */) {
      BOOST_FOREACH(s_elements_type const& elms, average_elements)
        add_scalar_obs(m, elms.get<0>(), is_signed);
      if (local_elements.size()) {
        BOOST_FOREACH(s_elements_type const& elms, local_elements)
          add_vector_obs(m, elms.get<0>(), site_label, is_signed);
      }
      if (correlation_elements.size()) {
        BOOST_FOREACH(p_elements_type const& elms, correlation_elements)
          add_vector_obs(m, elms.get<0>(), distance_label, is_signed);
      }
      if (strfactor_elements.size()) {
        BOOST_FOREACH(p_elements_type const& elms, strfactor_elements)
          add_vector_obs(m, elms.get<0>(), momenta_label, is_signed);
      }
    }

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
          // average
          BOOST_FOREACH(s_elements_type const& elms, emt.average_elements) {
            double v = 0;
            BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs,
              sites(lat.rg())) {
              int st = 0;
              BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs,
                sites(lat, rs))
                st += (1 - spins[vs]);
              v += elms.get<1>()[get(site_type_t(), lat.rg(), rs)][st];
            }
            m[elms.get<0>()] << sign * v / lat.volume();
          }

          // local
          std::valarray<double>& local = emt.local;
          BOOST_FOREACH(s_elements_type const& elms, emt.local_elements) {
            local = 0;
            BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs,
              sites(lat.rg())) {
              int st = 0;
              BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs,
                sites(lat, rs))
                st += (1 - spins[vs]);
              local[rs] = sign * elms.get<1>()[get(site_type_t(), lat.rg(), rs)][st];
            }
            m[elms.get<0>()] << local;
          }

          // correlation
          std::valarray<double>& corr = emt.corr;
          BOOST_FOREACH(p_elements_type const& elms, emt.correlation_elements) {
            corr = 0;
            if (emt.origin) {
              int origin = emt.origin.get();
              int st0 = 0;
              BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs0,
                sites(lat, origin))
                st0 += (1 - spins[vs0]);
              double v0 = elms.get<1>()[get(site_type_t(), lat.rg(), origin)][st0];
              BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs1,
                sites(lat.rg())) {
                int st1 = 0;
                BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs1,
                  sites(lat, rs1))
                  st1 += (1 - spins[vs1]);
                double v1 = elms.get<2>()[get(site_type_t(), lat.rg(), rs1)][st1];
                corr[rs1] = sign * v0 * v1;
              }
            } else {
              BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs0,
                sites(lat.rg())) {
                int st0 = 0;
                BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs0,
                  sites(lat, rs0))
                  st0 += (1 - spins[vs0]);
                double v0 = elms.get<1>()[get(site_type_t(), lat.rg(), rs0)][st0];
                BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs1,
                  sites(lat.rg())) {
                  int st1 = 0;
                  BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs1,
                    sites(lat, rs1))
                    st1 += (1 - spins[vs1]);
                  double v1 = elms.get<2>()[get(site_type_t(), lat.rg(), rs1)][st1];
                  int d = distance(lat, rs0, rs1);
                  corr[d] += v0 * v1;
                }
              }
              corr *= sign * emt.mltplcty;
            }
            m[elms.get<0>()] << corr;
          }

          // structure factor
          std::valarray<double>& sfac = emt.sfac;
          BOOST_FOREACH(p_elements_type const& elms, emt.strfactor_elements) {
            sfac = 0;
            int k = 0;
            typename momentum_iterator<lattice_t>::type mit, mit_end;
            for (boost::tie(mit, mit_end) = momenta(lat); mit != mit_end; ++mit, ++k) {
              std::complex<double> v0, v1;
              BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs,
                sites(lat.rg())) {
                int t = get(site_type_t(), lat.rg(), rs);
                int st = 0;
                BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs,
                  sites(lat, rs))
                  st += (1 - spins[vs]);
                v0 += elms.get<1>()[t][st] * mit.phase(emt.coordinate(rs));
                v1 += elms.get<2>()[t][st] * mit.phase(emt.coordinate(rs));
              }
              sfac[k] = std::real(std::conj(v0) * v1);
            }
            sfac *= (sign / lat.volume());
            m[elms.get<0>()] << sfac;
          }
        }
      };
    };
  };
};

} // end namespace looper

#endif // LOOPER_CUSTOM_H
