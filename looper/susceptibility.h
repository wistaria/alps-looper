/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include "measurement.h"
#include "divide_if_positive.h"

namespace looper {

struct susceptibility : public has_improved_estimator_tag, public has_normal_estimator_tag,
  public has_evaluator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<gauge_t, const typename lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool bipartite;
    gauge_map_t gauge;

    void initialize(alps::Parameters const&, lattice_t const& lat, bool /* is_signed */) {
      bipartite = is_bipartite(lat);
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      add_scalar_obs(m, "Magnetization", is_signed);
      add_scalar_obs(m, "Magnetization Density", is_signed);
      add_scalar_obs(m, "|Magnetization|", is_signed);
      add_scalar_obs(m, "|Magnetization Density|", is_signed);
      add_scalar_obs(m, "Magnetization^2", is_signed);
      add_scalar_obs(m, "Magnetization Density^2", is_signed);
      add_scalar_obs(m, "Magnetization^4", is_signed);
      add_scalar_obs(m, "Magnetization Density^4", is_signed);
      add_scalar_obs(m, "Susceptibility", is_signed);
      if (enable_improved_estimator) {
        add_scalar_obs(m, "Generalized Magnetization^2", is_signed);
        add_scalar_obs(m, "Generalized Magnetization Density^2", is_signed);
        add_scalar_obs(m, "Generalized Magnetization^4", is_signed);
        add_scalar_obs(m, "Generalized Magnetization Density^4", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility", is_signed);
      }
      if (bipartite) {
        add_scalar_obs(m, "Staggered Magnetization", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density", is_signed);
        add_scalar_obs(m, "|Staggered Magnetization|", is_signed);
        add_scalar_obs(m, "|Staggered Magnetization Density|", is_signed);
        add_scalar_obs(m, "Staggered Magnetization^2", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density^2", is_signed);
        add_scalar_obs(m, "Staggered Magnetization^4", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density^4", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility", is_signed);
        if (enable_improved_estimator) {
          add_scalar_obs(m, "Generalized Staggered Magnetization^2", is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization Density^2", is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization^4", is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization Density^4", is_signed);
          add_scalar_obs(m, "Generalized Staggered Susceptibility", is_signed);
        }
      }
    }
    
    struct improved_estimator {
      struct estimate {
        double usize0, umag0, usize, umag;
        double ssize0, smag0, ssize, smag;
        estimate() : usize0(0), umag0(0), usize(0), umag(0), ssize0(0), smag0(0), ssize(0),
          smag(0) {}
        void reset() {
          usize0 = umag0 = usize = umag = 0;
          ssize0 = smag0 = ssize = smag = 0;
        }
        estimate& operator+=(estimate const& rhs) {
          usize0 += rhs.usize0;
          umag0 += rhs.umag0;
          usize += rhs.usize;
          umag += rhs.umag;
          ssize0 += rhs.ssize0;
          smag0 += rhs.smag0;
          ssize += rhs.ssize;
          smag += rhs.smag;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -t, s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, double t, int s, int c) {
          usize += t * 0.5;
          umag  += t * (0.5-c);
          double gg = emt.gauge[s];
          ssize += gg * t * 0.5;
          smag  += gg * t * (0.5-c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
          usize0 += 0.5;
          umag0  += (0.5-c);
          double gg = emt.gauge[s];
          ssize0 += gg * 0.5;
          smag0  += gg * (0.5-c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
      };

      struct collector {
        double umag0, usize2, umag2, usize4, umag4, usize, umag;
        double smag0, ssize2, smag2, ssize4, smag4, ssize, smag;
        void reset() {
          umag0 = usize2 = umag2 = usize4 = umag4 = usize = umag = 0;
          smag0 = ssize2 = smag2 = ssize4 = smag4 = ssize = smag = 0;
        }
        collector& operator+=(collector const& coll) {
          umag0 += coll.umag0;
          usize2 += coll.usize2;
          umag2  += coll.umag2;
          usize4 += coll.usize4;
          umag4  += coll.umag4;
          usize  += coll.usize;
          umag   += coll.umag;
          smag0 += coll.smag0;
          ssize2 += coll.ssize2;
          smag2  += coll.smag2;
          ssize4 += coll.ssize4;
          smag4  += coll.smag4;
          ssize  += coll.ssize;
          smag   += coll.smag;
          return *this;
        }
        collector& operator+=(estimate const& est) {
          umag0  += est.umag0;
          usize2 += power2(est.usize0);
          umag2  += power2(est.umag0);
          usize4 += power4(est.usize0);
          umag4  += power4(est.umag0);
          usize  += power2(est.usize);
          umag   += power2(est.umag);
          smag0  += est.smag0;
          ssize2 += power2(est.ssize0);
          smag2  += power2(est.smag0);
          ssize4 += power4(est.ssize0);
          smag4  += power4(est.smag0);
          ssize  += power2(est.ssize);
          smag   += power2(est.smag);
          return *this;
        }
        template<typename M>
        void commit(M& m, lattice_t const& lat, double beta, double sign, int nop) const {
          double vol = lat.volume();
          m["Magnetization"] << 0.0;
          m["Magnetization Density"] << 0.0;
          m["|Magnetization|"] << sign * std::abs(umag0);
          m["|Magnetization Density|"] << sign * std::abs(umag0) / vol;
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
            m["|Staggered Magnetization|"] << sign * std::abs(smag0);
            m["|Staggered Magnetization Density|"] << sign * std::abs(smag0) / vol;
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
    };

    struct normal_estimator {
      struct collector {
        double umag, smag, umag_a, smag_a;
        collector() : umag(0), smag(0), umag_a(0), smag_a(0) {}
        void reset() { umag = smag = umag_a = smag_a = 0; }
        collector& operator+=(collector const& rhs) {
          umag += rhs.umag;
          smag += rhs.smag;
          umag_a += rhs.umag_a;
          smag_a += rhs.smag_a;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -t, s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, double t, int s, int c) {
          umag  += t * (0.5-c);
          smag  += emt.gauge[s] * t * (0.5-c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
          umag += 0.5-c;
          smag += (0.5-c) * emt.gauge[s];
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<typename M>
        void commit(M& m, lattice_t const& lat, double beta, double sign, int nop) const {
          double vol = lat.volume();
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
          m["Susceptibility"]
            << (typename is_sse<mc_type>::type() ?
                sign * beta * (dip(power2(umag_a), nop) + power2(umag)) / (nop + 1) / vol :
                sign * beta * power2(umag_a) / vol);
          if (is_bipartite(lat)) {
            m["Staggered Susceptibility"] 
              << (typename is_sse<mc_type>::type() ?
                  sign * beta * (dip(power2(smag_a), nop) + power2(smag)) / (nop + 1) / vol :
                  sign * beta * power2(smag_a) / vol);
          }
        }
      };
    };
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const&,
      alps::ObservableSet const& m_in) {
      try {
        alps::RealObsevaluator obse_m2 = m_in["Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Magnetization^4"];
        alps::RealObsevaluator eval("Binder Ratio of Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_m2 = m_in["Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Staggered Magnetization^4"];
        alps::RealObsevaluator eval("Binder Ratio of Staggered Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_m2 = m_in["Generalized Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Generalized Magnetization^4"];
        alps::RealObsevaluator eval("Binder Ratio of Generalized Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_m2 =
          m_in["Generalized Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 =
          m_in["Generalized Staggered Magnetization^4"];
        alps::RealObsevaluator eval("Binder Ratio of Generalized Staggered Magnetization");
        eval = power2(obse_m2) / obse_m4;
        m.addObservable(eval);
      } catch (...) {}
    }
  };
};

} // end namespace looper

#endif // LOOPER_SUSCEPTIBILITY_H
