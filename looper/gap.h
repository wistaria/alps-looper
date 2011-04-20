/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006-2011 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_GAP_H
#define LOOPER_GAP_H

#ifndef LOOPER_ONLY_PATH_INTEGRAL
# define LOOPER_ONLY_PATH_INTEGRAL
#endif

#include "measurement.h"
#include "correlation_length.h"
#include "time.h"
#include "type.h"
#ifdef HAVE_PARAPACK_13
# include <alps/math.hpp>
#else
# include <alps/numeric/is_nonzero.hpp>
#endif
#include <complex>

#ifdef HAVE_PARAPACK_13
using alps::is_nonzero;
#else
using alps::numeric::is_nonzero;
#endif

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct gap : public has_improved_estimator_tag, public has_normal_estimator_tag,
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

    void initialize(alps::Parameters const&, lattice_t const& lat, bool /* is_signed */,
      bool /* enable_improved_estimator */) {
      bipartite = is_bipartite(lat);
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      add_scalar_obs(m, "Susceptibility [w=2pi/beta]", is_signed);
      if (bipartite)
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]", is_signed);
      if (enable_improved_estimator) {
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]", is_signed);
        if (bipartite)
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=2pi/beta]", is_signed);
      }
    }

    struct improved_estimator {
      struct estimate {
        std::complex<double> upsize, upmag, spsize, spmag;
        estimate() : upsize(std::complex<double>(0,0)), upmag(std::complex<double>(0,0)),
          spsize(std::complex<double>(0,0)), spmag(std::complex<double>(0,0)) {
        }
        void reset(estimator_t const&) {
          upsize = upmag = spsize = spmag = std::complex<double>(0,0);
        }
        estimate& operator+=(estimate const& rhs) {
          upsize += rhs.upsize;
          upmag += rhs.upmag;
          spsize += rhs.spsize;
          spmag += rhs.spmag;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -ctime(t), s, c);
        }
        template<bool T>
        void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, -ctime(t), s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, ctime(t), s, c);
        }
        template<bool T>
        void end_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, ctime(t), s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, std::complex<double> const& ct, int s,
          int c) {
          upsize += ct * 0.5;
          upmag += ct * (0.5-c);
          double gg = emt.gauge[s];
          spsize += gg * ct * 0.5;
          spmag  += gg * ct * (0.5-c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start_bottom(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop_top(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
      };

      struct collector {
        double upsize2, upmag2, spsize2, spmag2;
        void reset(estimator_t const&) {
          upsize2 = 0;
          upmag2 = 0;
          spsize2 = 0;
          spmag2 = 0;
        }
        collector& operator+=(collector const& coll) {
          upsize2 += coll.upsize2;
          upmag2 += coll.upmag2;
          spsize2 += coll.spsize2;
          spmag2 += coll.spmag2;
          return *this;
        }
        collector& operator+=(estimate const& est) {
          upsize2 += power2(est.upsize);
          upmag2 += power2(est.upmag);
          spsize2 += power2(est.spsize);
          spmag2 += power2(est.spmag);
          return *this;
        }
        template<typename M>
        void commit(M& m, estimator_t const&, lattice_t const& lat, double beta, double sign,
          double) const {
          m["Susceptibility [w=2pi/beta]"] << sign * beta * upmag2 / power2(2*M_PI) / lat.volume();
          m["Generalized Susceptibility [w=2pi/beta]"]
            << sign * beta * upsize2 / power2(2*M_PI) / lat.volume();
          if (is_bipartite(lat)) {
            m["Staggered Susceptibility [w=2pi/beta]"]
              << sign * beta * spmag2 / power2(2*M_PI) / lat.volume();
            m["Generalized Staggered Susceptibility [w=2pi/beta]"]
              << sign * beta * spsize2 / power2(2*M_PI) / lat.volume();
          }
        }
      };
    };

    struct normal_estimator {
      struct collector {
        double umag, smag;
        std::complex<double> umag_a, smag_a;
        collector() : umag(0), smag(0), umag_a(std::complex<double>(0,0)),
          smag_a(std::complex<double>(0,0)) {}
        void reset(estimator_t const&) {
          umag = smag = 0;
          umag_a = smag_a = std::complex<double>(0,0);
        }
        collector& operator+=(collector const& rhs) {
          umag += rhs.umag;
          smag += rhs.smag;
          umag_a += rhs.umag_a;
          smag_a += rhs.smag_a;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -ctime(t), s, c);
        }
        template<bool T>
        void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, -ctime(t), s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, ctime(t), s, c);
        }
        template<bool T>
        void end_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, ctime(t), s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, std::complex<double> const& ct, int s,
          int c) {
          umag_a += ct * 0.5;
          smag_a += ct * (0.5-c) * emt.gauge[s];
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start_bottom(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop_top(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<typename M>
        void commit(M& m, estimator_t const&, lattice_t const& lat, double beta, double sign,
          double) const {
          m["Susceptibility [w=2pi/beta]"] <<
            sign * beta * power2(umag_a) / power2(2*M_PI) / lat.volume();
          if (is_bipartite(lat))
            m["Staggered Susceptibility [w=2pi/beta]"] <<
              sign * beta * power2(smag_a) / power2(2*M_PI) / lat.volume();
        }
      };
    };
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      double beta = 0;
      try {
        alps::RealObsevaluator beta_eval(m_in["Inverse Temperature"]);
        beta = beta_eval.mean();
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap [k=0]");
        alps::RealObsevaluator eval1("Gap [k=0]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap [k=pi]");
        alps::RealObsevaluator eval1("Gap [k=pi]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Gap");
        alps::RealObsevaluator eval1("Generalized Gap");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 =
          m_in["Generalized Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Staggered Gap");
        alps::RealObsevaluator eval1("Generalized Staggered Gap");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
    }
  };
};

struct gap4 : public has_improved_estimator_tag, public has_normal_estimator_tag,
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

    void initialize(alps::Parameters const&, lattice_t const& lat, bool /* is_signed */,
      bool /* enable_improved_estimator */) {
      bipartite = is_bipartite(lat);
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      add_scalar_obs(m, "Susceptibility [w=2pi/beta]", is_signed);
      add_scalar_obs(m, "Susceptibility [w=4pi/beta]", is_signed);
      if (bipartite)
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility [w=4pi/beta]", is_signed);
      if (enable_improved_estimator) {
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility [w=4pi/beta]", is_signed);
        if (bipartite) {
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=2pi/beta]", is_signed);
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=4pi/beta]", is_signed);
        }
      }
    }

    struct improved_estimator {
      struct estimate {
        std::complex<double> upsize, upmag, spsize, spmag;
        std::complex<double> up2size, up2mag, sp2size, sp2mag;
        estimate() : upsize(std::complex<double>(0,0)), upmag(std::complex<double>(0,0)),
          spsize(std::complex<double>(0,0)), spmag(std::complex<double>(0,0)),
          up2size(std::complex<double>(0,0)), up2mag(std::complex<double>(0,0)),
          sp2size(std::complex<double>(0,0)), sp2mag(std::complex<double>(0,0)) {
        }
        void reset(estimator_t const&) {
          upsize = upmag = spsize = spmag = std::complex<double>(0,0);
          up2size = up2mag = sp2size = sp2mag = std::complex<double>(0,0);
        }
        estimate& operator+=(estimate const& rhs) {
          upsize += rhs.upsize;
          upmag += rhs.upmag;
          spsize += rhs.spsize;
          spmag += rhs.spmag;
          up2size += rhs.up2size;
          up2mag += rhs.up2mag;
          sp2size += rhs.sp2size;
          sp2mag += rhs.sp2mag;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -ctime(t), -ctime2(t), s, c);
        }
        template<bool T>
        void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, -ctime(t), -ctime2(t), s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, ctime(t), ctime2(t), s, c);
        }
        template<bool T>
        void end_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, ctime(t), ctime2(t), s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, std::complex<double> const& ct,
          std::complex<double> const& ct2, int s, int c) {
          double gg = emt.gauge[s];
          upsize += ct * 0.5;
          upmag += ct * (0.5-c);
          spsize += gg * ct * 0.5;
          spmag  += gg * ct * (0.5-c);
          up2size += ct2 * 0.5;
          up2mag += ct2 * (0.5-c);
          sp2size += gg * ct2 * 0.5;
          sp2mag  += gg * ct2 * (0.5-c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start_bottom(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop_top(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
    };

      struct collector {
        double upsize2, upmag2, spsize2, spmag2;
        double up2size2, up2mag2, sp2size2, sp2mag2;
        void reset(estimator_t const&) {
          upsize2 = 0;
          upmag2 = 0;
          spsize2 = 0;
          spmag2 = 0;
          up2size2 = 0;
          up2mag2 = 0;
          sp2size2 = 0;
          sp2mag2 = 0;
        }
        collector& operator+=(collector const& coll) {
          upsize2 += coll.upsize2;
          upmag2 += coll.upmag2;
          spsize2 += coll.spsize2;
          spmag2 += coll.spmag2;
          up2size2 += coll.up2size2;
          up2mag2 += coll.up2mag2;
          sp2size2 += coll.sp2size2;
          sp2mag2 += coll.sp2mag2;
          return *this;
        }
        collector& operator+=(estimate const& est) {
          upsize2 += power2(est.upsize);
          upmag2 += power2(est.upmag);
          spsize2 += power2(est.spsize);
          spmag2 += power2(est.spmag);
          up2size2 += power2(est.up2size);
          up2mag2 += power2(est.up2mag);
          sp2size2 += power2(est.sp2size);
          sp2mag2 += power2(est.sp2mag);
          return *this;
        }
        template<typename M>
        void commit(M& m, estimator_t const&, lattice_t const& lat, double beta, double sign,
          double) const {
          m["Susceptibility [w=2pi/beta]"]
            << sign * beta * upmag2 / power2(2*M_PI) / lat.volume();
          m["Susceptibility [w=4pi/beta]"]
            << sign * beta * up2mag2 / power2(4*M_PI) / lat.volume();
          m["Generalized Susceptibility [w=2pi/beta]"]
            << sign * beta * upsize2 / power2(2*M_PI) / lat.volume();
          m["Generalized Susceptibility [w=4pi/beta]"]
            << sign * beta * up2size2 / power2(4*M_PI) / lat.volume();
          if (is_bipartite(lat)) {
            m["Staggered Susceptibility [w=2pi/beta]"]
              << sign * beta * spmag2 / power2(2*M_PI) / lat.volume();
            m["Staggered Susceptibility [w=4pi/beta]"]
              << sign * beta * sp2mag2 / power2(4*M_PI) / lat.volume();
            m["Generalized Staggered Susceptibility [w=2pi/beta]"]
              << sign * beta * spsize2 / power2(2*M_PI) / lat.volume();
            m["Generalized Staggered Susceptibility [w=4pi/beta]"]
              << sign * beta * sp2size2 / power2(4*M_PI) / lat.volume();
          }
        }
      };
    };

    struct normal_estimator {
      struct collector {
        double umag, smag;
        std::complex<double> umag_a, smag_a;
        std::complex<double> u2mag_a, s2mag_a;
        collector() : umag(0), smag(0), umag_a(std::complex<double>(0,0)),
          smag_a(std::complex<double>(0,0)), u2mag_a(std::complex<double>(0,0)),
          s2mag_a(std::complex<double>(0,0)) {}
        void reset(estimator_t const&) {
          umag = smag = 0;
          umag_a = smag_a = std::complex<double>(0,0);
          u2mag_a = s2mag_a = std::complex<double>(0,0);
        }
        collector& operator+=(collector const& rhs) {
          umag += rhs.umag;
          smag += rhs.smag;
          umag_a += rhs.umag_a;
          smag_a += rhs.smag_a;
          u2mag_a += rhs.u2mag_a;
          s2mag_a += rhs.s2mag_a;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -ctime(t), -ctime2(t), s, c);
        }
        template<bool T>
        void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, -ctime(t), -ctime2(t), s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, ctime(t), ctime2(t), s, c);
        }
        template<bool T>
        void end_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, ctime(t), ctime2(t), s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, std::complex<double> const& ct,
          std::complex<double> const& ct2, int s, int c) {
          double gg = emt.gauge[s];
          umag_a += ct * 0.5;
          smag_a += gg * ct * (0.5-c);
          u2mag_a += ct2 * 0.5;
          s2mag_a += gg * ct2 * (0.5-c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start_bottom(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop_top(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }

        void end_s(estimator_t const& emt, lattice_t const&, std::complex<double> const& ct, int s,
          int c) {
        }
        template<typename M>
        void commit(M& m, estimator_t const&, lattice_t const& lat, double beta, double sign,
          double) const {
          m["Susceptibility [w=2pi/beta]"] <<
            sign * beta * power2(umag_a) / power2(2*M_PI) / lat.volume();
          m["Susceptibility [w=4pi/beta]"] <<
            sign * beta * power2(u2mag_a) / power2(4*M_PI) / lat.volume();
          if (is_bipartite(lat)) {
            m["Staggered Susceptibility [w=2pi/beta]"] <<
              sign * beta * power2(smag_a) / power2(2*M_PI) / lat.volume();
            m["Staggered Susceptibility [w=4pi/beta]"] <<
              sign * beta * power2(s2mag_a) / power2(4*M_PI) / lat.volume();
          }
        }
      };
    };
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      double beta = 0;
      try {
        alps::RealObsevaluator beta_eval(m_in["Inverse Temperature"]);
        beta = beta_eval.mean();
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap [k=0]");
        alps::RealObsevaluator eval1("Gap [k=0]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator obse_s4 = m_in["Susceptibility [w=4pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap (4th moment estimator) [k=0]");
        alps::RealObsevaluator eval1("Gap (4th moment estimator) [k=0]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
          eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap [k=pi]");
        alps::RealObsevaluator eval1("Gap [k=pi]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator obse_s4 = m_in["Staggered Susceptibility [w=4pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap (4th moment estimator) [k=pi]");
        alps::RealObsevaluator eval1("Gap (4th moment estimator) [k=pi]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
          eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Gap");
        alps::RealObsevaluator eval1("Generalized Gap");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator obse_s4 = m_in["Generalized Susceptibility [w=4pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Gap (4th moment estimator)");
        alps::RealObsevaluator eval1("Generalized Gap (4th moment estimator)");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
          eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 =
          m_in["Generalized Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Staggered Gap");
        alps::RealObsevaluator eval1("Generalized Staggered Gap");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 =
          m_in["Generalized Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator obse_s4 =
          m_in["Generalized Staggered Susceptibility [w=4pi/beta]"];
        alps::RealObsevaluator
          eval0("Inverse Generalized Staggered Gap (4th moment estimator)");
        alps::RealObsevaluator
          eval1("Generalized Staggered Gap (4th moment estimator)");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
          eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
    }
  };
};

struct generalized_gap4_and_correlation_length : public has_improved_estimator_tag,
  public has_normal_estimator_tag, public has_evaluator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<gauge_t, const typename lattice_t::virtual_graph_type,
      double>::type gauge_map_t;

    double dq_abs;
    gauge_map_t gauge;
    std::vector<std::complex<double> > phase1;

    void initialize(alps::Parameters const& params, lattice_t const& lat, bool /* is_signed */,
      bool /* enable_improved_estimator */) {
      if (!is_bipartite(lat)) {
        std::cerr << "Error: gap4_and_correlation_length: lattice is not bipartite\n";
        boost::throw_exception(std::invalid_argument("gap4_and_correlation_length: lattice is not bipartite"));
      }
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);

      int dim = lat.graph_helper().dimension();
      int nvs = num_sites(lat.vg());

      std::vector<double> dq(dim, 0.0);

      if (!params.defined("DQ_OVER_TWO_PI")) {
        std::cerr << "parameter DQ_OVER_TWO_PI (dq/2pi) is not defined\n";
        boost::throw_exception(std::invalid_argument("correlation_length"));
      }

      dq = parse_vec(static_cast<std::string>(params["DQ_OVER_TWO_PI"]), params);
      if (dq.size() != dim) {
        std::cerr << "inconsistent dimension of DQ_OVER_TWO_PI (dq/2pi)\n";
        boost::throw_exception(std::invalid_argument("correlation_length::initialize"));
      }

      phase1.resize(nvs);
      typename alps::graph_helper<typename lattice_t::real_graph_type>::vector_type coord;
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg())) {
        coord = lat.graph_helper().coordinate(get(real_site_t(), lat.vg(), s));
        double p1 = 0;
        for (int i = 0; i < dim; ++i) {
          p1 += dq[i] * coord[i];
        }
        phase1[s] = std::complex<double>(std::cos(p1), std::sin(p1));
      }

      dq_abs = 0;
      for (int i = 0; i < dim; ++i) dq_abs += dq[i] * dq[i];
      dq_abs = std::sqrt(dq_abs);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      if (enable_improved_estimator) {
        add_scalar_obs(m, "Generalized Susceptibility", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility [w=4pi/beta]", is_signed);
        add_scalar_obs(m, "Generalized Spin Dynamic Structure Factor at (q,w) = (dq,0)",
          is_signed);
        add_scalar_obs(m, "|dq|");
      } else {
        add_scalar_obs(m, "Staggered Susceptibility", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility [w=4pi/beta]", is_signed);
        add_scalar_obs(m, "Staggered Spin Dynamic Structure Factor at (q,w) = (dq,0)",
          is_signed);
      }
      add_scalar_obs(m, "|dq|");
    }

    struct improved_estimator {
      struct estimate {
        std::complex<double> upsize, up2size, sq1;
        double usize;
        estimate() : upsize(std::complex<double>(0,0)), up2size(std::complex<double>(0,0)),
          sq1(std::complex<double>(0,0)), usize(0) {
        }
        void reset(estimator_t const&) {
          upsize = up2size = sq1 = std::complex<double>(0,0);
          usize = 0;
        }
        estimate& operator+=(estimate const& rhs) {
          upsize += rhs.upsize;
          up2size += rhs.up2size;
          sq1 += rhs.sq1;
          usize += rhs.usize;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -t, -ctime(t), -ctime2(t), s, c);
        }
        template<bool T>
        void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, -t, -ctime(t), -ctime2(t), s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, ctime(t), ctime2(t), s, c);
        }
        template<bool T>
        void end_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, ctime(t), ctime2(t), s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, double t,
          std::complex<double> const& ct, std::complex<double> const& ct2, int s, int c) {
          upsize += ct;
          up2size += ct2;
          sq1 += emt.phase1[s] * t;
          usize += t;
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start_bottom(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop_top(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
      };

      struct collector {
        double upsize2, up2size2, str1, usize2;
        void reset(estimator_t const&) {
          upsize2 = 0;
          up2size2 = 0;
          str1 = 0;
          usize2 = 0;
        }
        collector& operator+=(collector const& coll) {
          upsize2 += coll.upsize2;
          up2size2 += coll.up2size2;
          str1 += coll.str1;
          usize2 += coll.usize2;
          return *this;
        }
        collector& operator+=(estimate const& est) {
          upsize2 += power2(est.upsize);
          up2size2 += power2(est.up2size);
          str1 += power2(est.sq1);
          usize2 += power2(est.usize);
          return *this;
        }
        template<typename M>
        void commit(M& m, estimator_t const& emt, lattice_t const& lat, double beta, double sign,
          double) const {
          m["|dq|"] << emt.dq_abs;
          m["Generalized Susceptibility [w=2pi/beta]"]
            << sign * 0.25 * beta * upsize2 / power2(2*M_PI) / lat.volume();
          m["Generalized Susceptibility [w=4pi/beta]"]
            << sign * 0.25 * beta * up2size2 / power2(4*M_PI) / lat.volume();
          m["Generalized Spin Dynamic Structure Factor at (q,w) = (dq,0)"]
            << sign * 0.25 * beta * str1 / lat.volume();
          m["Generalized Susceptibility"]
            << sign * 0.25 * beta * usize2 / lat.volume();
#ifdef STD_OUTPUT
        std::cout << sign * 0.25 * beta * upsize2 / power2(2*M_PI) / lat.volume() << ' '
                  << sign * 0.25 * beta * up2size2 / power2(4*M_PI) / lat.volume() << ' '
                  << sign * 0.25 * beta * str1 / lat.volume() << ' '
                  << sign * 0.25 * beta * usize2 / lat.volume() << ' ';
#endif
        }
      };
    };
    struct normal_estimator {
      struct collector {
        double smag_a;
        std::complex<double> qmag_a, q2mag_a, sq1;
        collector() : smag_a(0), qmag_a(std::complex<double>(0,0)),
          q2mag_a(std::complex<double>(0,0)), sq1(std::complex<double>(0,0)) {}
        void reset(estimator_t const&) {
          smag_a = 0;
          qmag_a = std::complex<double>(0,0);
          q2mag_a = std::complex<double>(0,0);
          sq1 = std::complex<double>(0,0);
        }
        collector& operator+=(collector const& rhs) {
          smag_a += rhs.smag_a;
          qmag_a += rhs.qmag_a;
          q2mag_a += rhs.q2mag_a;
          sq1 += rhs.sq1;
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, -t, -ctime(t), -ctime2(t), s, c);
        }
        template<bool T>
        void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, -t, -ctime(t), -ctime2(t), s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, ctime(t), ctime2(t), s, c);
        }
        template<bool T>
        void end_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, ctime(t), ctime2(t), s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const&, double t,
          std::complex<double> const& ct, std::complex<double> const& ct2, int s, int c) {
          double gg = emt.gauge[s];
          smag_a += gg * t * (0.5-c);
          qmag_a += gg * ct * (0.5-c);
          q2mag_a += gg * ct2 * (0.5-c);
          sq1 += gg * emt.phase1[s] * t * (0.5-c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void end_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start_bottom(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          begin_s(emt, lat, t, s, c);
        }
        template<bool T>
        void start(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          begin_s(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
          int c) {
          end_s(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<bool T>
        void stop_top(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t,
          int s, int c) {
          end_s(emt, lat, t, s, c);
        }
        template<typename M>
        void commit(M& m, estimator_t const& emt, lattice_t const& lat, double beta, double sign,
          double) const {
          m["|dq|"] << emt.dq_abs;
          m["Staggered Susceptibility [w=2pi/beta]"] <<
            sign * beta * power2(qmag_a) / power2(2*M_PI) / lat.volume();
          m["Staggered Susceptibility [w=4pi/beta]"] <<
            sign * beta * power2(q2mag_a) / power2(4*M_PI) / lat.volume();
          m["Staggered Spin Dynamic Structure Factor at (q,w) = (dq,0)"] <<
            sign * beta * power2(sq1) / lat.volume();
          m["Staggered Susceptibility"] <<
            sign * beta * power2(smag_a) / lat.volume();
#ifdef STD_OUTPUT
        std::cout
          << sign * beta * power2(qmag_a) / power2(2*M_PI) / lat.volume()
          << sign * beta * power2(q2mag_a) / power2(4*M_PI) / lat.volume()
          << sign * beta * power2(sq1) / lat.volume()
          << sign * beta * power2(smag_a) / lat.volume();
#endif
        }
      };
    };
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      double beta = 0;
      try {
        alps::RealObsevaluator beta_eval(m_in["Inverse Temperature"]);
        beta = beta_eval.mean();
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Gap");
        alps::RealObsevaluator eval1("Generalized Gap");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator obse_s4 = m_in["Generalized Susceptibility [w=4pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Generalized Gap (4th moment estimator)");
        alps::RealObsevaluator eval1("Generalized Gap (4th moment estimator)");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
          eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator dq_abs_eval(m_in["|dq|"]);
        alps::RealObsevaluator obse_s0(m_in["Generalized Susceptibility"]);
        alps::RealObsevaluator obse_s1(m_in["Generalized Spin Dynamic Structure Factor at (q,w) = (dq,0)"]);
        double dq_abs = dq_abs_eval.mean();
        alps::RealObsevaluator eval("Correlation Length");
        eval = sqrt(obse_s0 / obse_s1 - 1) / dq_abs;
        m.addObservable(eval);
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap [k=pi]");
        alps::RealObsevaluator eval1("Gap [k=pi]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean())) {
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
        alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
        alps::RealObsevaluator obse_s4 = m_in["Staggered Susceptibility [w=4pi/beta]"];
        alps::RealObsevaluator eval0("Inverse Gap (4th moment estimator) [k=pi]");
        alps::RealObsevaluator eval1("Gap (4th moment estimator) [k=pi]");
        if (alps::numeric::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
          eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      } catch (...) {}
      try {
        alps::RealObsevaluator dq_abs_eval(m_in["|dq|"]);
        alps::RealObsevaluator obse_s0(m_in["Staggered Susceptibility"]);
        alps::RealObsevaluator obse_s1(m_in["Staggered Spin Dynamic Structure Factor at (q,w) = (dq,0)"]);
        double dq_abs = dq_abs_eval.mean();
        alps::RealObsevaluator eval("Correlation Length");
        eval = sqrt(obse_s0 / obse_s1 - 1) / dq_abs;
        m.addObservable(eval);
      } catch (...) {}
    }
  };
};

} // end namespace looper

#endif // LOOPER_GAP_H
