/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
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
#include "time.h"
#include "type.h"
#include <alps/math.hpp>
#include <complex>

namespace looper {

struct gap : public has_evaluator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<gauge_t, const typename lattice_t::virtual_graph_type,
      double>::type gauge_map_t;

    bool bipartite, improved;
    gauge_map_t gauge;

    void initialize(alps::Parameters const& /* params */, lattice_t const& lat,
      bool /* is_signed */, bool use_improved_estimator) {
      bipartite = is_bipartite(lat);
      improved = use_improved_estimator;
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      add_scalar_obs(m, "Susceptibility [w=2pi/beta]", is_signed);
      if (bipartite)
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]", is_signed);
      if (improved) {
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]", is_signed);
        if (bipartite)
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=2pi/beta]", is_signed);
      }
    }

    // improved estimator

    struct estimate {
      std::complex<double> upsize, upmag, spsize, spmag;
      estimate() : upsize(std::complex<double>(0,0)), upmag(std::complex<double>(0,0)),
        spsize(std::complex<double>(0,0)), spmag(std::complex<double>(0,0)) {
      }
      void init() {
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
      void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
                   int c) {
        end_s(emt, lat, -ctime(t), s, c);
      }
      void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      template<bool T>
      void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
                    int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      template<bool T>
      void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
                    int s, int c) {
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
    void init_estimate(estimate& est) const { est.init(); }

    struct collector {
      double upsize2, upmag2, spsize2, spmag2;
      void init() {
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
      void commit(M& m, lattice_t const& lat, double beta, int, double sign) const {
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
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& /* spins */, std::vector<OP> const& operators,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators, std::vector<int>& spins_c) {
      if (!typename is_path_integral<mc_type>::type()) return;
      if (improved) return;

      double umag = 0;
      double smag = 0;
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg())) {
        umag += (0.5-spins[s]);
        smag += (0.5-spins[s]) * gauge[s];
      }
      std::complex<double> umag_a(0, 0);
      std::complex<double> smag_a(0, 0);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          const std::complex<double> p = ctime(oi->time());
          umag_a += p * umag;
          smag_a += p * smag;
          if (oi->is_site()) {
            const unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += (1-2*spins_c[s]);
            smag += gauge[s] * (1-2*spins_c[s]);
          } else {
            const unsigned int s0 = source(oi->pos(), lat.vg());
            const unsigned int s1 = target(oi->pos(), lat.vg());
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            umag += (1-2*spins_c[s0]) + (1-2*spins_c[s1]);
            smag += gauge[s0] * (1-2*spins_c[s0]) + gauge[s1] * (1-2*spins_c[s1]);
          }
          umag_a -= p * umag;
          smag_a -= p * smag;
        }
      }
      m["Susceptibility [w=2pi/beta]"] <<
        sign * beta * power2(umag_a) / power2(2*M_PI) / lat.volume();
      if (is_bipartite(lat))
        m["Staggered Susceptibility [w=2pi/beta]"] <<
          sign * beta * power2(smag_a) / power2(2*M_PI) / lat.volume();
    }
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      if (m_in.has("Inverse Temperature")) {
        alps::RealObsevaluator beta_eval(m_in["Inverse Temperature"]);
        if (beta_eval.count() == 0) return;
        double beta = beta_eval.mean();
        if (m_in.has("Susceptibility") &&
            m_in.has("Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Gap [k=0]");
            alps::RealObsevaluator eval1("Gap [k=0]");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Staggered Susceptibility") &&
            m_in.has("Staggered Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Gap [k=pi]");
            alps::RealObsevaluator eval1("Gap [k=pi]");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Generalized Susceptibility") &&
            m_in.has("Generalized Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Generalized Gap");
            alps::RealObsevaluator eval1("Generalized Gap");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Generalized Staggered Susceptibility") &&
            m_in.has("Generalized Staggered Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 =
            m_in["Generalized Staggered Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Generalized Staggered Gap");
            alps::RealObsevaluator eval1("Generalized Staggered Gap");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
      }
    }
  };
};

struct gap4 : public has_evaluator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<gauge_t, const typename lattice_t::virtual_graph_type,
      double>::type gauge_map_t;

    bool bipartite, improved;
    gauge_map_t gauge;

    void initialize(alps::Parameters const& /* params */, lattice_t const& lat,
      bool /* is_signed */, bool use_improved_estimator) {
      bipartite = is_bipartite(lat);
      improved = use_improved_estimator;
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      add_scalar_obs(m, "Susceptibility [w=2pi/beta]", is_signed);
      add_scalar_obs(m, "Susceptibility [w=4pi/beta]", is_signed);
      if (bipartite) {
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility [w=4pi/beta]", is_signed);
      }
      if (improved) {
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility [w=4pi/beta]", is_signed);
        if (bipartite) {
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=2pi/beta]", is_signed);
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=4pi/beta]", is_signed);
        }
      }
    }

    // improved estimator

    struct estimate {
      std::complex<double> upsize, upmag, spsize, spmag;
      std::complex<double> up2size, up2mag, sp2size, sp2mag;
      estimate() : upsize(std::complex<double>(0,0)), upmag(std::complex<double>(0,0)),
                   spsize(std::complex<double>(0,0)), spmag(std::complex<double>(0,0)),
                   up2size(std::complex<double>(0,0)), up2mag(std::complex<double>(0,0)),
                   sp2size(std::complex<double>(0,0)), sp2mag(std::complex<double>(0,0)) {
      }
      void init() {
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
      void begin_s(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int s,
                   int c) {
        end_s(emt, lat, -ctime(t), -ctime2(t),s, c);
      }
      void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      template<bool T>
      void begin_bs(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
                    int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      template<bool T>
      void begin_bt(estimator_t const& emt, lattice_t const& lat, imaginary_time<T> const& t, int,
                    int s, int c) {
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
        const double gg = emt.gauge[s];
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
    void init_estimate(estimate& est) const { est.init(); }

    struct collector {
      double upsize2, upmag2, spsize2, spmag2;
      double up2size2, up2mag2, sp2size2, sp2mag2;
      void init() {
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
      void commit(M& m, lattice_t const& lat, double beta, int, double sign) const {
        m["Susceptibility [w=2pi/beta]"] << sign * beta * upmag2 / power2(2*M_PI) / lat.volume();
        m["Susceptibility [w=4pi/beta]"] << sign * beta * up2mag2 / power2(4*M_PI) / lat.volume();
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
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& /* spins */, std::vector<OP> const& operators,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators, std::vector<int>& spins_c) {
      if (!typename is_path_integral<mc_type>::type()) return;
      if (improved) return;

      double umag = 0;
      double smag = 0;
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type s, sites(lat.vg())) {
        umag += (0.5-spins[s]);
        smag += (0.5-spins[s]) * gauge[s];
      }
      std::complex<double> umag_a(0, 0);
      std::complex<double> smag_a(0, 0);
      std::complex<double> u2mag_a(0, 0);
      std::complex<double> s2mag_a(0, 0);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          const std::complex<double> p = ctime(oi->time());
          const std::complex<double> p2 = p * p;
          umag_a += p * umag;
          smag_a += p * smag;
          u2mag_a += p2 * umag;
          s2mag_a += p2 * smag;
          if (oi->is_site()) {
            const unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += (1-2*spins_c[s]);
            smag += gauge[s] * (1-2*spins_c[s]);
          } else {
            const unsigned int s0 = source(oi->pos(), lat.vg());
            const unsigned int s1 = target(oi->pos(), lat.vg());
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            umag += (1-2*spins_c[s0]) + (1-2*spins_c[s1]);
            smag += gauge[s0] * (1-2*spins_c[s0]) + gauge[s1] * (1-2*spins_c[s1]);
          }
          umag_a -= p * umag;
          smag_a -= p * smag;
          u2mag_a -= p2 * umag;
          s2mag_a -= p2 * smag;
        }
      }
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

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      if (m_in.has("Inverse Temperature")) {
        alps::RealObsevaluator beta_eval(m_in["Inverse Temperature"]);
        if (beta_eval.count() == 0) return;
        double beta = beta_eval.mean();
        if (m_in.has("Susceptibility") &&
            m_in.has("Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Gap [k=0]");
            alps::RealObsevaluator eval1("Gap [k=0]");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Susceptibility") &&
            m_in.has("Susceptibility [w=2pi/beta]") &&
            m_in.has("Susceptibility [w=4pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Susceptibility [w=2pi/beta]"];
          alps::RealObsevaluator obse_s4 = m_in["Susceptibility [w=4pi/beta]"];
          if (obse_s0.count() && obse_s2.count() &&
              alps::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
            alps::RealObsevaluator eval0("Inverse Gap (4th moment estimator) [k=0]");
            alps::RealObsevaluator eval1("Gap (4th moment estimator) [k=0]");
            eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Staggered Susceptibility") &&
            m_in.has("Staggered Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Gap [k=pi]");
            alps::RealObsevaluator eval1("Gap [k=pi]");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Staggered Susceptibility") &&
            m_in.has("Staggered Susceptibility [w=2pi/beta]") &&
            m_in.has("Staggered Susceptibility [w=4pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Staggered Susceptibility [w=2pi/beta]"];
          alps::RealObsevaluator obse_s4 = m_in["Staggered Susceptibility [w=4pi/beta]"];
          if (obse_s0.count() && obse_s2.count() &&
              alps::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
            alps::RealObsevaluator eval0("Inverse Gap (4th moment estimator) [k=pi]");
            alps::RealObsevaluator eval1("Gap (4th moment estimator) [k=pi]");
            eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Generalized Susceptibility") &&
            m_in.has("Generalized Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Generalized Gap");
            alps::RealObsevaluator eval1("Generalized Gap");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Generalized Susceptibility") &&
            m_in.has("Generalized Susceptibility [w=2pi/beta]") &&
            m_in.has("Generalized Susceptibility [w=4pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
          alps::RealObsevaluator obse_s2 = m_in["Generalized Susceptibility [w=2pi/beta]"];
          alps::RealObsevaluator obse_s4 = m_in["Generalized Susceptibility [w=4pi/beta]"];
          if (obse_s0.count() && obse_s2.count() &&
              alps::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
            alps::RealObsevaluator eval0("Inverse Generalized Gap (4th moment estimator)");
            alps::RealObsevaluator eval1("Generalized Gap (4th moment estimator)");
            eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Generalized Staggered Susceptibility") &&
            m_in.has("Generalized Staggered Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 =
            m_in["Generalized Staggered Susceptibility [w=2pi/beta]"];
          if (obse_s0.count() && obse_s2.count() && alps::is_nonzero<1>(obse_s2.mean())) {
            alps::RealObsevaluator eval0("Inverse Generalized Staggered Gap");
            alps::RealObsevaluator eval1("Generalized Staggered Gap");
            eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
        if (m_in.has("Generalized Staggered Susceptibility") &&
            m_in.has("Generalized Staggered Susceptibility [w=2pi/beta]") &&
            m_in.has("Generalized Staggered Susceptibility [w=4pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 =
            m_in["Generalized Staggered Susceptibility [w=2pi/beta]"];
          alps::RealObsevaluator obse_s4 =
            m_in["Generalized Staggered Susceptibility [w=4pi/beta]"];
          if (obse_s0.count() && obse_s2.count() &&
              alps::is_nonzero<1>(obse_s2.mean() - obse_s4.mean())) {
            alps::RealObsevaluator
              eval0("Inverse Generalized Staggered Gap (4th moment estimator)");
            alps::RealObsevaluator
              eval1("Generalized Staggered Gap (4th moment estimator)");
            eval0 = sqrt(3.0*(obse_s0-obse_s2)/(obse_s2-obse_s4) - 1) / (4*M_PI/beta);
            eval1 = 1.0 / eval0;
            m.addObservable(eval0);
            m.addObservable(eval1);
          }
        }
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_GAP_H
