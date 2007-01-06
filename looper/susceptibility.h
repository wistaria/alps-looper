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

#ifndef LOOPER_SUSCEPTIBILITY_H
#define LOOPER_SUSCEPTIBILITY_H

#include "divide_if_positive.h"
#include "measurement.h"
#include "type.h"

namespace looper {

struct susceptibility {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<gauge_t,
              const typename lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool improved;
    gauge_map_t gauge;

    template<typename M>
    void initialize(M& m, alps::Parameters const& /* params */, lattice_t const& lat,
      bool is_signed, bool use_improved_estimator) {
      improved = use_improved_estimator;
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);

      add_scalar_obs(m, "Magnetization", is_signed);
      add_scalar_obs(m, "Magnetization Density", is_signed);
      add_scalar_obs(m, "|Magnetization|", is_signed);
      add_scalar_obs(m, "|Magnetization Density|", is_signed);
      add_scalar_obs(m, "Magnetization^2", is_signed);
      add_scalar_obs(m, "Magnetization Density^2", is_signed);
      add_scalar_obs(m, "Magnetization^4", is_signed);
      add_scalar_obs(m, "Magnetization Density^4", is_signed);
      add_scalar_obs(m, "Susceptibility", is_signed);
      if (use_improved_estimator) {
        add_scalar_obs(m, "Generalized Magnetization^2", is_signed);
        add_scalar_obs(m, "Generalized Magnetization Density^2", is_signed);
        add_scalar_obs(m, "Generalized Magnetization^4", is_signed);
        add_scalar_obs(m, "Generalized Magnetization Density^4", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility", is_signed);
      }
      if (is_bipartite(lat)) {
        add_scalar_obs(m, "Staggered Magnetization", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density", is_signed);
        add_scalar_obs(m, "|Staggered Magnetization|", is_signed);
        add_scalar_obs(m, "|Staggered Magnetization Density|", is_signed);
        add_scalar_obs(m, "Staggered Magnetization^2", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density^2", is_signed);
        add_scalar_obs(m, "Staggered Magnetization^4", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density^4", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility", is_signed);
        if (use_improved_estimator) {
          add_scalar_obs(m, "Generalized Staggered Magnetization^2",
                         is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization Density^2",
                         is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization^4",
                         is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization Density^4",
                         is_signed);
          add_scalar_obs(m, "Generalized Staggered Susceptibility",
                         is_signed);
        }
      }
    }

    // improved estimator

    struct estimate {
      gauge_map_t gauge;
      double usize0, umag0, usize, umag;
      double ssize0, smag0, ssize, smag;
      void init(gauge_map_t map) {
        gauge = map;
        usize0 = 0;
        umag0 = 0;
        usize = 0;
        umag = 0;
        ssize0 = 0;
        smag0 = 0;
        ssize = 0;
        smag = 0;
      }
      void start_s(lattice_t const& lat, double t, int s, int c) { term_s(lat, -t, s, c); }
      void start_bs(lattice_t const& lat, double t, int, int s, int c) { start_s(lat, t, s, c); }
      void start_bt(lattice_t const& lat, double t, int, int s, int c) { start_s(lat, t, s, c); }
      void term_s(lattice_t const&, double t, int s, int c) {
        usize += t * 0.5;
        umag  += t * (0.5-c);
        double gg = gauge[s];
        ssize += gg * t * 0.5;
        smag  += gg * t * (0.5-c);
      }
      void term_bs(lattice_t const& lat, double t, int, int s, int c) { term_s(lat, t, s, c); }
      void term_bt(lattice_t const& lat, double t, int, int s, int c) { term_s(lat, t, s, c); }
      void at_bot(lattice_t const& lat, double t, int s, int c) {
        start_s(lat, t, s, c);
        usize0 += 0.5;
        umag0  += (0.5-c);
        double gg = gauge[s];
        ssize0 += gg * 0.5;
        smag0  += gg * (0.5-c);
      }
      void at_top(lattice_t const& lat, double t, int s, int c) { term_s(lat, t, s, c); }
    };
    void init_estimate(estimate& est) const { est.init(gauge); }

    struct collector {
      double usize2, umag2, usize4, umag4, usize, umag;
      double ssize2, smag2, ssize4, smag4, ssize, smag;
      void init() {
        usize2 = 0; umag2 = 0; usize4 = 0; umag4 = 0;
        usize = 0; umag = 0;
        ssize2 = 0; smag2 = 0; ssize4 = 0; smag4 = 0;
        ssize = 0; smag = 0;
      }
      template<typename EST>
      collector operator+(EST const& est) {
        usize2 += power2(est.usize0);
        umag2  += power2(est.umag0);
        usize4 += power4(est.usize0);
        umag4  += power4(est.umag0);
        usize  += power2(est.usize);
        umag   += power2(est.umag);
        ssize2 += power2(est.ssize0);
        smag2  += power2(est.smag0);
        ssize4 += power4(est.ssize0);
        smag4  += power4(est.smag0);
        ssize  += power2(est.ssize);
        smag   += power2(est.smag);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const& lat, double beta, int nop, double sign) const {
        double vol = lat.volume();
        m["Magnetization"] << 0.0;
        m["Magnetization Density"] << 0.0;
        m["|Magnetization|"] << 0.0;
        m["|Magnetization Density|"] << 0.0;
        m["Magnetization^2"] << sign * umag2;
        m["Magnetization Density^2"] << sign * umag2 / power2(vol);
        m["Magnetization^4"] << sign * (3 * power2(umag2) - 2 * umag4);
        m["Magnetization Density^4"]
          << sign * (3 * power2(umag2) - 2 * umag4) / power4(vol);
        m["Susceptibility"]
          << (typename is_sse<mc_type>::type() ?
              sign * beta * (dip(umag, nop) + umag2) / (nop + 1) / vol :
              sign * beta * umag / vol);
        m["Generalized Magnetization^2"] << sign * usize2;
        m["Generalized Magnetization Density^2"]
          << sign * usize2 / power2(vol);
        m["Generalized Magnetization^4"]
          << sign * (3 * power2(usize2) - 2 * usize4);
        m["Generalized Magnetization Density^4"]
          << sign * (3 * power2(usize2) - 2 * usize4) / power4(vol);
        m["Generalized Susceptibility"]
          << (typename is_sse<mc_type>::type() ?
              sign * beta * (dip(usize, nop) + usize2) / (nop + 1) / vol :
              sign * beta * usize / vol);
        if (is_bipartite(lat)) {
          m["Staggered Magnetization"] << 0.0;
          m["Staggered Magnetization Density"] << 0.0;
          m["|Staggered Magnetization|"] << 0.0;
          m["|Staggered Magnetization Density|"] << 0.0;
          m["Staggered Magnetization^2"] << sign * smag2;
          m["Staggered Magnetization Density^2"] << sign * smag2 / power2(vol);
          m["Staggered Magnetization^4"]
            << sign * (3 * power2(smag2) - 2 * smag4);
          m["Staggered Magnetization Density^4"]
            << sign * (3 * power2(smag2) - 2 * smag4) / power4(vol);
          m["Staggered Susceptibility"]
            << (typename is_sse<mc_type>::type() ?
                sign * beta * (dip(smag, nop) + smag2) / (nop + 1) / vol :
                sign * beta * smag / vol);
          m["Generalized Staggered Magnetization^2"] << sign * ssize2;
          m["Generalized Staggered Magnetization Density^2"]
            << sign * ssize2 / power2(vol);
          m["Generalized Staggered Magnetization^4"]
            << sign * (3 * power2(ssize2) - 2 * ssize4);
          m["Generalized Staggered Magnetization Density^4"]
            << sign * (3 * power2(ssize2) - 2 * ssize4) / power4(vol);
          m["Generalized Staggered Susceptibility"]
            << (typename is_sse<mc_type>::type() ?
                sign * beta * (dip(ssize, nop) + ssize2) / (nop + 1) / vol :
                sign * beta * ssize / vol);
        }
      }
    };
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& /* spins */, std::vector<OP> const& operators,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators,
      std::vector<int>& spins_c) {
      if (improved) return;

      double vol = lat.volume();
      double nop = operators.size();
      double umag = 0;
      double smag = 0;
      typename virtual_site_iterator<lattice_t>::type si, si_end;
      for (boost::tie(si, si_end) = sites(lat.vg()); si != si_end; ++si) {
        umag += 0.5-spins[*si];
        smag += (0.5-spins[*si]) * gauge[*si];
      }
      m["Magnetization"] << sign * umag;
      m["Magnetization Density"] << sign * umag / vol;
      m["|Magnetization|"] << sign * std::abs(umag);
      m["|Magnetization Density|"] << sign * std::abs(umag) / vol;
      m["Magnetization^2"] << sign * power2(umag);
      m["Magnetization Density^2"] << sign * power2(umag / vol);
      m["Magnetization^4"] << sign * power4(umag);
      m["Magnetization Density^4"] << sign * power4(umag / vol);
      if (is_bipartite(lat)) {
        m["Staggered Magnetization"] << sign * smag;
        m["Staggered Magnetization Density"] << sign * smag / vol;
        m["|Staggered Magnetization|"] << sign * std::abs(smag);
        m["|Staggered Magnetization Density|"] << sign * std::abs(smag) / vol;
        m["Staggered Magnetization^2"] << sign * power2(smag);
        m["Staggered Magnetization Density^2"] << sign * power2(smag / vol);
        m["Staggered Magnetization^4"] << sign * power4(smag);
        m["Staggered Magnetization Density^4"] << sign * power4(smag / vol);
      }
      double umag_a = 0; /* 0 * umag; */
      double smag_a = 0; /* 0 * smag; */
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      double t = 0;
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          proceed(typename is_path_integral<mc_type>::type(), t, *oi);
          umag_a += t * umag;
          smag_a += t * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += 1-2*spins_c[s];
            smag += gauge[s] * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = source(oi->pos(), lat.vg());
            unsigned int s1 = target(oi->pos(), lat.vg());
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            umag += 1-2*spins_c[s0] + 1-2*spins_c[s1];
            smag += gauge[s0] * (1-2*spins_c[s0])
              + gauge[s1] * (1-2*spins_c[s1]);
          }
          umag_a -= t * umag;
          smag_a -= t * smag;
        }
        proceed(typename is_sse<mc_type>::type(), t);
      }
      if (typename is_path_integral<mc_type>::type()) {
        umag_a += umag;
        m["Susceptibility"] << sign * beta * power2(umag_a) / vol;
        smag_a += smag;
        if (is_bipartite(lat))
          m["Staggered Susceptibility"] << sign * beta * power2(smag_a) / vol;
      } else {
        umag_a += nop * umag;
        m["Susceptibility"]
          << sign * beta * (dip(power2(umag_a), nop) + power2(umag))
          / (nop + 1) / vol;
        smag_a += nop * smag;
        if (is_bipartite(lat))
          m["Staggered Susceptibility"]
            << sign * beta * (dip(power2(smag_a), nop) + power2(smag))
            / (nop + 1) / vol;
      }
    }
  };

  struct evaluator
  {
    static void evaluate(alps::ObservableSet& m,
                         alps::Parameters const& /* params */,
                         alps::ObservableSet const& m_in)
    {
      if (m_in.has("Magnetization^2") && m_in.has("Magnetization^4")) {
        alps::RealObsevaluator obse_m2 = m_in["Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Magnetization^4"];
        alps::RealObsevaluator eval("Binder Ratio of Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      }
      if (m_in.has("Staggered Magnetization^2") &&
          m_in.has("Staggered Magnetization^4")) {
        alps::RealObsevaluator obse_m2 = m_in["Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Staggered Magnetization^4"];
        alps::RealObsevaluator eval("Binder Ratio of Staggered Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      }
      if (m_in.has("Generalized Magnetization^2") &&
          m_in.has("Generalized Magnetization^4")) {
        alps::RealObsevaluator obse_m2 = m_in["Generalized Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Generalized Magnetization^4"];
        alps::RealObsevaluator
          eval("Binder Ratio of Generalized Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      }
      if (m_in.has("Generalized Staggered Magnetization^2") &&
          m_in.has("Generalized Staggered Magnetization^4")) {
        alps::RealObsevaluator obse_m2 =
          m_in["Generalized Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 =
          m_in["Generalized Staggered Magnetization^4"];
        alps::RealObsevaluator
          eval("Binder Ratio of Generalized Staggered Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_SUSCEPTIBILITY_H
