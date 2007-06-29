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

#ifndef LOOPER_LOCALSUS_MEASUREMENT_H
#define LOOPER_LOCALSUS_MEASUREMENT_H

#include "measurement.h"
#include "type.h"

namespace looper {

struct local_susceptibility {
  typedef dumb_measurement<local_susceptibility> dumb;

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<real_site_t,
              const typename lattice_t::virtual_graph_type,
              typename real_site_descriptor<lattice_t>::type>::type
              real_site_map_t;
    typedef typename alps::property_map<gauge_t,
              const typename lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool measure;
    bool measure_local;
    bool measure_type;
    bool improved;
    real_site_map_t real_site;
    gauge_map_t gauge;

    unsigned int types;
    std::valarray<double> type_nrs;

    mutable int next_id;

    mutable std::vector<int> vs2c; // table of virtual site to cluster
    mutable std::vector<double> usize;
    mutable std::vector<double> umag;
    mutable std::vector<double> ssize;
    mutable std::vector<double> smag;

    std::valarray<double> tlumag, tlumag_a, tlsmag, tlsmag_a; // temporary
    std::valarray<double> tumag, tumag_a, tsmag, tsmag_a; // temporary

    template<typename M>
    void initialize(M& m, alps::Parameters const& params, lattice_t const& lat,
      bool is_signed, bool use_improved_estimator) {
      measure_local = params.value_or_default("MEASURE[Local Susceptibility]", false);
      measure_type = params.value_or_default("MEASURE[Site Type Susceptibility]", false);
      measure = measure_local || measure_type;

      improved = use_improved_estimator;
      real_site = alps::get_or_default(real_site_t(), lat.vg(), 0);
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);

      if (measure_type) {
        types = 0;
        BOOST_FOREACH(int s, sites(lat.rg())) {
          int t = get(site_type_t(), lat.rg(), s);
          if (t < 0)
            boost::throw_exception(std::runtime_error("negative site type"));
          if (t >= types) types = t + 1;
        }
        type_nrs.resize(types); type_nrs = 0;
        BOOST_FOREACH(int s, sites(lat.rg())) type_nrs[get(site_type_t(), lat.rg(), s)] += 1;
      }

      next_id = 0;

      vs2c.resize(num_sites(lat.vg()));
      usize.resize(num_sites(lat.vg()));
      umag.resize(num_sites(lat.vg()));
      ssize.resize(num_sites(lat.vg()));
      smag.resize(num_sites(lat.vg()));

      if (measure_local) {
        tlumag.resize(num_sites(lat.rg()));
        tlumag_a.resize(num_sites(lat.rg()));
        tlsmag.resize(num_sites(lat.rg()));
        tlsmag_a.resize(num_sites(lat.rg()));

        add_vector_obs(m, "Local Magnetization", is_signed);
        add_vector_obs(m, "Local Susceptibility", is_signed);
        add_vector_obs(m, "Local Field Susceptibility", is_signed);
        if (is_bipartite(lat))
          add_vector_obs(m, "Staggered Local Susceptibility", is_signed);
        if (improved) {
          add_vector_obs(m, "Generalized Local Susceptibility", is_signed);
          add_vector_obs(m, "Generalized Local Field Susceptibility", is_signed);
          if (is_bipartite(lat))
            add_vector_obs(m, "Generalized Staggered Local Susceptibility", is_signed);
        }
      }
      if (measure_type) {
        tumag.resize(type_nrs.size());
        tumag_a.resize(type_nrs.size());
        tsmag.resize(type_nrs.size());
        tsmag_a.resize(type_nrs.size());

        add_vector_obs(m, "Number of Sites of Each Type", is_signed);
        add_vector_obs(m, "Site Type Magnetization", is_signed);
        add_vector_obs(m, "Site Type Magnetization Density", is_signed);
        add_vector_obs(m, "Site Type Susceptibility", is_signed);
        add_vector_obs(m, "Site Type Field Susceptibility", is_signed);
        if (use_improved_estimator) {
          add_vector_obs(m, "Generalized Site Type Susceptibility", is_signed);
        }
        if (is_bipartite(lat)) {
          add_vector_obs(m, "Staggered Site Type Magnetization", is_signed);
          add_vector_obs(m, "Staggered Site Type Magnetization Density", is_signed);
          add_vector_obs(m, "Staggered Site Type Susceptibility", is_signed);
          add_vector_obs(m, "Staggered Site Type Field Susceptibility", is_signed);
          if (use_improved_estimator) {
            add_vector_obs(m, "Generalized Staggered Site Type Susceptibility", is_signed);
          }
        }
      }
    }

    // improved estimator

    struct estimate {
      bool cross_itb;
      int cluster_id;
      gauge_map_t gauge;
      int* id_ptr;
      std::vector<int>* vs2c_ptr;
      std::vector<double>* usize_ptr;
      std::vector<double>* umag_ptr;
      std::vector<double>* ssize_ptr;
      std::vector<double>* smag_ptr;
      void init(gauge_map_t map, int* id, std::vector<int>* vs2c,
        std::vector<double>* usize, std::vector<double>* umag,
        std::vector<double>* ssize, std::vector<double>* smag) {
        cross_itb = false;
        cluster_id = -1;
        gauge = map;
        id_ptr = id;
        vs2c_ptr = vs2c;
        usize_ptr = usize;
        umag_ptr = umag;
        ssize_ptr = ssize;
        smag_ptr = smag;
      }
      void start_s(lattice_t const& lat, double t, int s, int c) { term_s(lat, -t, s, c); }
      void start_bs(lattice_t const& lat, double t, int, int s, int c) { start_s(lat, t, s, c); }
      void start_bt(lattice_t const& lat, double t, int, int s, int c) { start_s(lat, t, s, c); }
      void term_s(lattice_t const&, double t, int s, int c) {
        if (cross_itb) {
          (*usize_ptr)[cluster_id] += t * 0.5;
          (*umag_ptr)[cluster_id] += t * (0.5 - c);
          (*ssize_ptr)[cluster_id] += gauge[s] * t * 0.5;
          (*smag_ptr)[cluster_id] += gauge[s] * t * (0.5 - c);
        }
      }
      void term_bs(lattice_t const& lat, double t, int, int s, int c) { term_s(lat, t, s, c); }
      void term_bt(lattice_t const& lat, double t, int, int s, int c) { term_s(lat, t, s, c); }
      void at_bot(lattice_t const& lat, double t, int s, int c) {
        if (!cross_itb) {
          cross_itb = true;
          cluster_id = *id_ptr;
          (*id_ptr) += 1;
        }
        (*vs2c_ptr)[s] = cluster_id;
        (*usize_ptr)[cluster_id] = 0;
        (*umag_ptr)[cluster_id] = 0;
        (*ssize_ptr)[cluster_id] = 0;
        (*smag_ptr)[cluster_id] = 0;
        start_s(lat, t, s, c);
      }
      void at_top(lattice_t const& lat, double t, int s, int c) { term_s(lat, t, s, c); }
    };
    void init_estimate(estimate& est) const {
      est.init(gauge, &next_id, &vs2c, &usize, &umag, &ssize, &smag);
    }

    typedef typename dumb::template estimator<MC, LAT, TIME>::collector collector;
    void init_collector(collector const&) const {}

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& /* operators */,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& /* coll */) {

      if (!typename looper::is_path_integral<mc_type>::type()) return;
      if (!measure) return;
      if (!improved) return;

      int nrs = num_vertices(lat.rg());
      next_id = 0;

      if (measure_local) {
        m["Local Magnetization"] << std::valarray<double>(0., nrs);

        tlumag = 0;
        tlumag_a = 0;
        tlsmag = 0;
        tlsmag_a = 0;
        BOOST_FOREACH(int s, sites(lat.vg())) {
          int r = real_site[s];
          tlumag[r] += umag[vs2c[s]] * (0.5 - spins[s]);
          tlumag_a[r] += usize[vs2c[s]] * 0.5;
          tlsmag[r] += gauge[s] * smag[vs2c[s]] * (0.5 - spins[s]);
          tlsmag_a[r] += gauge[s] * ssize[vs2c[s]] * 0.5;
        }
        m["Local Susceptibility"] << std::valarray<double>(sign * beta * tlumag);
        m["Generalized Local Susceptibility"] << std::valarray<double>(sign * beta * tlumag_a);
        if (is_bipartite(lat)) {
          m["Staggered Local Susceptibility"] << std::valarray<double>(sign * beta * tlsmag);
          m["Generalized Staggered Local Susceptibility"]
            << std::valarray<double>(sign * beta * tlsmag_a);
        }
      }

      if (measure_type) {
        m["Number of Sites of Each Type"] << type_nrs;

        tumag = 0;
        tumag_a = 0;
        tsmag = 0;
        tsmag_a = 0;
        BOOST_FOREACH(int s, sites(lat.vg())) {
          int p = get(site_type_t(), lat.rg(), real_site[s]);
          tumag[p] += umag[vs2c[s]] * (0.5 - spins[s]);
          tumag_a[p] += usize[vs2c[s]] * 0.5;
          tsmag[p] += gauge[s] * smag[vs2c[s]] * (0.5 - spins[s]);
          tsmag_a[p] += gauge[s] * ssize[vs2c[s]] * 0.5;
        }
        m["Site Type Susceptibility"] << std::valarray<double>(sign * beta * tumag / type_nrs);
        m["Generalized Site Type Susceptibility"]
          << std::valarray<double>(sign * beta * tumag_a / type_nrs);
        if (is_bipartite(lat)) {
          m["Staggered Site Type Susceptibility"]
            << std::valarray<double>(sign * beta * tsmag / type_nrs);
          m["Generalized Staggered Site Type Susceptibility"]
            << std::valarray<double>(sign * beta * tsmag_a / type_nrs);
        }
      }
    }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators,
      std::vector<int>& spins_c) {
      if (!typename looper::is_path_integral<mc_type>::type()) return;
      if (!measure) return;

      int nrs = num_vertices(lat.rg());

      if (measure_local) {
        tlumag = 0;
        tlsmag = 0;
        double umag = 0;
        double smag = 0;
        BOOST_FOREACH(int s, sites(lat.vg())) {
          int r = real_site[s];
          tlumag[r] += 0.5-spins[s];
          tlsmag[r] += gauge(s) * (0.5-spins[s]);
          umag += 0.5-spins[s];
          smag += gauge[s] * (0.5-spins[s]);
        }

        tlumag_a = 0; /* 0 * tlumag; */
        tlsmag_a = 0; /* 0 * tlsmag; */
        double umag_a = 0; /* 0 * umag; */
        double smag_a = 0; /* 0 * smag; */
        std::copy(spins.begin(), spins.end(), spins_c.begin());
        double t = 0;
        BOOST_FOREACH(OP const& op, operators) {
          if (op.is_offdiagonal()) {
            proceed(typename is_path_integral<mc_type>::type(), t, op);
            if (op.is_site()) {
              int s = op.pos();
              int r = real_site[s];
              tlumag_a[r] += t * tlumag[r];
              tlsmag_a[r] += t * tlsmag[r];
              umag_a += t * umag;
              smag_a += t * smag;
              spins_c[s] ^= 1;
              tlumag[r] += 1-2*spins_c[s];
              tlsmag[r] += gauge[s] * (1-2*spins_c[s]);
              umag += 1-2*spins_c[s];
              smag += gauge[s] * (1-2*spins_c[s]);
              tlumag_a[r] -= t * tlumag[r];
              tlsmag_a[r] -= t * tlsmag[r];
              umag_a -= t * umag;
              smag_a -= t * smag;
            } else {
              int s0 = source(op.pos(), lat.vg());
              int s1 = target(op.pos(), lat.vg());
              int r0 = real_site[s0];
              int r1 = real_site[s1];
              tlumag_a[r0] += t * tlumag[r0];
              tlumag_a[r1] += t * tlumag[r1];
              tlsmag_a[r0] += t * tlsmag[r0];
              tlsmag_a[r1] += t * tlsmag[r1];
              umag_a += t * umag;
              smag_a += t * smag;
              spins_c[s0] ^= 1;
              spins_c[s1] ^= 1;
              tlumag[r0] += 1-2*spins_c[s0];
              tlumag[r1] += 1-2*spins_c[s1];
              tlsmag[r0] += gauge[s0] * (1-2*spins_c[s0]);
              tlsmag[r1] += gauge[s1] * (1-2*spins_c[s1]);
              umag += (1-2*spins_c[s0]) + (1-2*spins_c[s1]);
              smag += (gauge[s0] * (1-2*spins_c[s0])) + (gauge[s1] * (1-2*spins_c[s1]));
              tlumag_a[r0] -= t * tlumag[r0];
              tlumag_a[r1] -= t * tlumag[r1];
              tlsmag_a[r0] -= t * tlsmag[r0];
              tlsmag_a[r1] -= t * tlsmag[r1];
              umag_a -= t * umag;
              smag_a -= t * smag;
            }
          }
          proceed(typename is_sse<mc_type>::type(), t);
        }
        for (int r = 0; r < nrs; ++r) {
          tlumag_a[r] += tlumag[r];
          tlsmag_a[r] += tlsmag[r];
        }
        umag_a += umag;
        smag_a += smag;

        m["Local Magnetization"] << std::valarray<double>(sign * tlumag_a);
        m["Local Field Susceptibility"]
          << std::valarray<double>(sign * beta * power2(tlumag_a));

        if (!improved) {
          m["Local Susceptibility"]
            << std::valarray<double>(sign * beta * umag_a * tlumag_a);
          if (is_bipartite(lat))
            m["Staggered Local Susceptibility"]
              << std::valarray<double>(sign * beta * smag_a * tlsmag_a);
        }
      }

      if (measure_type) {
        tumag = 0;
        tsmag = 0;
        double umag = 0;
        double smag = 0;
        BOOST_FOREACH(int s, sites(lat.vg())) {
          int p = get(site_type_t(), lat.rg(), real_site[s]);
          tumag[p] += 0.5 - spins[s];
          tsmag[p] += gauge[s] * (0.5 - spins[s]);
          umag += 0.5 - spins[s];
          smag += gauge[s] * (0.5 - spins[s]);
        }

        tumag_a = 0;
        tsmag_a = 0;
        double umag_a = 0;
        double smag_a = 0;
        std::copy(spins.begin(), spins.end(), spins_c.begin());
        double t = 0;
        BOOST_FOREACH(OP const& op, operators) {
          if (op.is_offdiagonal()) {
            proceed(typename is_path_integral<mc_type>::type(), t, op);
            if (op.is_site()) {
              int s = op.pos();
              int p = get(site_type_t(), lat.rg(), real_site[s]);
              tumag_a[p] += t * tumag[p];
              tsmag_a[p] += t * tsmag[p];
              umag_a += t * umag;
              smag_a += t * smag;
              spins_c[s] ^= 1;
              tumag[p] += 1-2*spins_c[s];
              tsmag[p] += gauge[s] * (1-2*spins_c[s]);
              umag += 1-2*spins_c[s];
              smag += gauge[s] * (1-2*spins_c[s]);
              tumag_a[p] -= t * tumag[p];
              tsmag_a[p] -= t * tsmag[p];
              umag_a -= t * umag;
              smag_a -= t * smag;
            } else {
              int s0 = source(op.pos(), lat.vg());
              int s1 = target(op.pos(), lat.vg());
              int p0 = get(site_type_t(), lat.rg(), real_site[s0]);
              int p1 = get(site_type_t(), lat.rg(), real_site[s1]);
              if (p0 == p1) {
                tumag_a[p0] += t * tumag[p0];
                tsmag_a[p0] += t * tsmag[p0];
              } else {
                tumag_a[p0] += t * tumag[p0];
                tumag_a[p1] += t * tumag[p1];
                tsmag_a[p0] += t * tsmag[p0];
                tsmag_a[p1] += t * tsmag[p1];
              }
              umag_a += t * umag;
              smag_a += t * smag;
              spins_c[s0] ^= 1;
              spins_c[s1] ^= 1;
              tumag[p0] += 1-2*spins_c[s0];
              tumag[p1] += 1-2*spins_c[s1];
              tsmag[p0] += gauge[s0] * (1-2*spins_c[s0]);
              tsmag[p1] += gauge[s1] * (1-2*spins_c[s1]);
              umag += (1-2*spins_c[s0]) + (1-2*spins_c[s1]);
              smag += (gauge[s0] * (1-2*spins_c[s0])) +
                (gauge[s1] * (1-2*spins_c[s1]));
              if (p0 == p1) {
                tumag_a[p0] -= t * tumag[p0];
                tsmag_a[p0] -= t * tsmag[p0];
              } else {
                tumag_a[p0] -= t * tumag[p0];
                tumag_a[p1] -= t * tumag[p1];
                tsmag_a[p0] -= t * tsmag[p0];
                tsmag_a[p1] -= t * tsmag[p1];
              }
              umag_a -= t * umag;
              smag_a -= t * smag;
            }
          }
          proceed(typename is_sse<mc_type>::type(), t);
        }
        for (int p = 0; p < types; ++p) {
          tumag_a[p] += tumag[p];
          tsmag_a[p] += tsmag[p];
        }
        umag_a += umag;
        smag_a += smag;

        m["Site Type Magnetization"] << std::valarray<double>(sign * tumag_a);
        m["Site Type Magnetization Density"] << std::valarray<double>(sign * tumag_a / type_nrs);
        m["Site Type Field Susceptibility"]
          << std::valarray<double>(sign * beta * power2(tumag_a) / type_nrs);
        if (is_bipartite(lat)) {
          m["Staggered Site Type Magnetization"] << std::valarray<double>(sign * tsmag_a);
          m["Staggered Site Type Magnetization Density"]
            << std::valarray<double>(sign * tsmag_a / type_nrs);
          m["Staggered Site Type Field Susceptibility"]
            << std::valarray<double>(sign * beta * power2(tsmag_a) / type_nrs);
        }
        if (!improved) {
          m["Number of Sites of Each Type"] << type_nrs;
          m["Site Type Susceptibility"]
            << std::valarray<double>(sign * beta * umag_a * tumag_a / type_nrs);
          if (is_bipartite(lat)) {
            m["Staggered Site Type Susceptibility"]
              << std::valarray<double>(sign * beta * smag_a * tsmag_a / type_nrs);
          }
        }
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_TYPEDSUS_MEASUREMENT_H
