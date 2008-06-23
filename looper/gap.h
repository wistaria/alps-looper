/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006-2008 by Synge Todo <wistaria@comp-phys.org>
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
#include <complex>

namespace looper {

struct gap : public has_evaluator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<gauge_t, const typename lattice_t::virtual_graph_type,
      double>::type gauge_map_t;

    bool improved;
    gauge_map_t gauge;

    template<typename M>
    void initialize(M& m, alps::Parameters const& /* params */, lattice_t const& lat,
      bool is_signed, bool use_improved_estimator) {
      improved = use_improved_estimator;
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
      add_scalar_obs(m, "Susceptibility [w=2pi/beta]", is_signed);
      if (is_bipartite(lat))
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]", is_signed);
      if (use_improved_estimator) {
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]", is_signed);
        if (is_bipartite(lat))
          add_scalar_obs(m, "Generalized Staggered Susceptibility [w=2pi/beta]", is_signed);
      }
    }

    // improved estimator

    struct estimate {
      gauge_map_t gauge;
      std::complex<double> upsize, upmag, spsize, spmag;
      void init(gauge_map_t map) {
        gauge = map;
        upsize = std::complex<double>(0,0);
        upmag = std::complex<double>(0,0);
        spsize = std::complex<double>(0,0);
        spmag = std::complex<double>(0,0);
      }
      void start_s(lattice_t const& lat, double t, int s, int c) { term_s(lat, -ctime(t), s, c); }
      template<typename T>
      void start_s(lattice_t const& lat, imaginary_time<T> const& t, int s, int c) {
        term_s(lat, -ctime(t), s, c);
      }
      void start_bs(lattice_t const& lat, double t, int, int s, int c) { start_s(lat, t, s, c); }
      void start_bt(lattice_t const& lat, double t, int, int s, int c) { start_s(lat, t, s, c); }
      template<typename T>
      void start_bs(lattice_t const& lat, imaginary_time<T> const& t, int, int s, int c) {
        start_s(lat, t, s, c);
      }
      template<typename T>
      void start_bt(lattice_t const& lat, imaginary_time<T> const& t, int, int s, int c) {
        start_s(lat, t, s, c);
      }
      void term_s(lattice_t const& lat, double t, int s, int c) { term_s(lat, ctime(t), s, c); }
      template<typename T>
      void term_s(lattice_t const& lat, imaginary_time<T> const& t, int s, int c) {
        term_s(lat, ctime(t), s, c);
      }
      void term_s(lattice_t const&, std::complex<double> const& ct, int s, int c) {
        upsize += ct * 0.5;
        upmag += ct * (0.5-c);
        double gg = gauge[s];
        spsize += gg * ct * 0.5;
        spmag  += gg * ct * (0.5-c);
      }
      void term_bs(lattice_t const& lat, double t, int, int s, int c) { term_s(lat, t, s, c); }
      void term_bt(lattice_t const& lat, double t, int, int s, int c) { term_s(lat, t, s, c); }
      template<typename T>
      void term_bs(lattice_t const& lat, imaginary_time<T> const& t, int, int s, int c) {
        term_s(lat, t, s, c);
      }
      template<typename T>
      void term_bt(lattice_t const& lat, imaginary_time<T> const& t, int, int s, int c) {
        term_s(lat, t, s, c);
      }
      void at_bot(lattice_t const& lat, double t, int s, int c) { start_s(lat, t, s, c); }
      template<typename T>
      void at_bot(lattice_t const& lat, imaginary_time<T> const& t, int s, int c) {
        start_s(lat, t, s, c);
      }
      void at_top(lattice_t const& lat, double t, int s, int c) { term_s(lat, t, s, c); }
      template<typename T>
      void at_top(lattice_t const& lat, imaginary_time<T> const& t, int s, int c) {
        term_s(lat, t, s, c);
      }
    };
    void init_estimate(estimate& est) const { est.init(gauge); }

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
          std::complex<double> p = ctime(oi->time());
          umag_a += p * umag;
          smag_a += p * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += (1-2*spins_c[s]);
            smag += gauge[s] * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = source(oi->pos(), lat.vg());
            unsigned int s1 = target(oi->pos(), lat.vg());
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
          if (obse_s0.count() && obse_s2.count()) {
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
          if (obse_s0.count() && obse_s2.count()) {
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
          if (obse_s0.count() && obse_s2.count()) {
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
          if (obse_s0.count() && obse_s2.count()) {
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

} // end namespace looper

#endif // LOOPER_GAP_H
