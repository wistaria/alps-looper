/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_GAP_MEASUREMENT_H
#define LOOPER_GAP_MEASUREMENT_H

#include "measurement.h"
#include "type.h"
#include <complex>

namespace looper {

struct gap
{
  template<typename MC, typename LAT, typename TIME>
  struct estimator
  {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<gauge_t,
              const typename lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    gauge_map_t gauge;

    template<typename M>
    void initialize(M& m, alps::Parameters const& /* params */,
                    lattice_t const& lat,
                    bool is_signed, bool use_improved_estimator)
    {
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
      if (is_bipartite(lat))
        add_scalar_obs(m, "Staggered Susceptibility [w=2pi/beta]",
                       is_signed);
      if (use_improved_estimator)
        add_scalar_obs(m, "Generalized Susceptibility [w=2pi/beta]",
                       is_signed);
    }

    // improved estimator

    struct estimate
    {
      std::complex<double> p;
      void init() { p = std::complex<double>(0,0); }
      void start_s(lattice_t const&, double t, int, int)
      { p -= ctime(t); }
      void start_s(lattice_t const&,
        imaginary_time<boost::mpl::true_> const& t, int, int)
      { p -= t.ctime_; }
      void start_bs(lattice_t const& lat, double t, int, int s, int c)
      { start_s(lat, t, s, c); }
      void start_bt(lattice_t const& lat, double t, int, int s, int c)
      { start_s(lat, t, s, c); }
      void start_bs(lattice_t const& lat,
        imaginary_time<boost::mpl::true_> const& t, int, int s , int c)
      { start_s(lat, t, s, c); }
      void start_bt(lattice_t const& lat,
        imaginary_time<boost::mpl::true_> const& t, int, int s , int c)
      { start_s(lat, t, s, c); }
      void term_s(lattice_t const&, double t, int, int)
      { p += ctime(t); }
      void term_s(lattice_t const&,
        imaginary_time<boost::mpl::true_> const& t, int, int)
      { p += t.ctime_; }
      void term_bs(lattice_t const& lat, double t, int, int s, int c)
      { term_s(lat, t, s, c); }
      void term_bt(lattice_t const& lat, double t, int, int s, int c)
      { term_s(lat, t, s, c); }
      void term_bs(lattice_t const& lat,
        imaginary_time<boost::mpl::true_> const& t, int, int s, int c)
      { term_s(lat, t, s, c); }
      void term_bt(lattice_t const& lat,
        imaginary_time<boost::mpl::true_> const& t, int, int s, int c)
      { term_s(lat, t, s, c); }
      void at_bot(lattice_t const& lat, double t, int s, int c)
      { start_s(lat, t, s, c); }
      void at_bot(lattice_t const& lat,
        imaginary_time<boost::mpl::true_> const& t, int s, int c)
      { start_s(lat, t, s, c); }
      void at_top(lattice_t const& lat, double t, int s, int c)
      { term_s(lat, t, s, c); }
      void at_top(lattice_t const& lat,
        imaginary_time<boost::mpl::true_> const& t, int s, int c)
      { term_s(lat, t, s, c); }
    };
    void init_estimate(estimate& est) const { est.init(); }

    struct collector
    {
      double p;
      void init() { p = 0; }
      template<typename EST>
      collector operator+(EST const& cm) {
        p += power2(cm.p);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const& lat, double beta, int,
                  double sign) const {
        m["Generalized Susceptibility [w=2pi/beta]"]
          << sign * beta * p / power2(4*M_PI) / num_sites(lat.rg());
      }
    };
    void init_collector(collector& coll) const {
      coll.init();
    }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m,
                              lattice_t const& lat,
                              double beta, double sign,
                              std::vector<int> const& /* spins */,
                              std::vector<OP> const& operators,
                              std::vector<int> const& /* spins_c */,
                              std::vector<FRAGMENT> const& /* fragments */,
                              collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat,
                            double beta, double sign,
                            std::vector<int> const& spins,
                            std::vector<OP> const& operators,
                            std::vector<int>& spins_c) {
      if (!typename is_path_integral<mc_type>::type()) return;
      if (!is_bipartite(lat)) return;

      int nrs = num_sites(lat.rg());

      double smag = 0;
      typename virtual_site_iterator<lattice_t>::type si, si_end;
      for (boost::tie(si, si_end) = sites(lat.vg()); si != si_end; ++si)
        smag += (0.5-spins[*si]) * gauge[*si];
      std::complex<double> smag_a(0, 0);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          std::complex<double> p = ctime(oi->time());
          smag_a += p * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            smag += gauge[s] * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = source(oi->pos(), lat.vg());
            unsigned int s1 = target(oi->pos(), lat.vg());
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            smag += gauge[s0] * (1-2*spins_c[s0])
              + gauge[s1] * (1-2*spins_c[s1]);
          }
          smag_a -= p * smag;
        }
      }
      m["Staggered Susceptibility [w=2pi/beta]"] <<
        sign * beta * power2(smag_a) / power2(2*M_PI) / nrs;
    }
  };

  struct evaluator
  {
    static void evaluate(alps::ObservableSet& m,
                         alps::Parameters const& /* params */,
                         alps::ObservableSet const& m_in)
    {
      if (m_in.has("Inverse Temperature")) {
        double beta =
          alps::RealObsevaluator(m_in["Inverse Temperature"]).mean();
        if (m_in.has("Generalized Susceptibility") &&
            m_in.has("Generalized Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Generalized Susceptibility"];
          alps::RealObsevaluator obse_s2 =
            m_in["Generalized Susceptibility [w=2pi/beta]"];
          alps::RealObsevaluator eval0("Inverse Gap");
          alps::RealObsevaluator eval1("Gap");
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
        if (m_in.has("Staggered Susceptibility") &&
            m_in.has("Staggered Susceptibility [w=2pi/beta]")) {
          alps::RealObsevaluator obse_s0 = m_in["Staggered Susceptibility"];
          alps::RealObsevaluator obse_s2 =
            m_in["Staggered Susceptibility [w=2pi/beta]"];
          alps::RealObsevaluator eval0("Inverse Gap [k=pi]");
          alps::RealObsevaluator eval1("Gap [k=pi]");
          eval0 = sqrt(obse_s0/obse_s2 - 1) / (2*M_PI/beta);
          eval1 = 1.0 / eval0;
          m.addObservable(eval0);
          m.addObservable(eval1);
        }
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_GAP_MEASUREMENT_H
