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

#ifndef GAP_MEASUREMENT_H
#define GAP_MEASUREMENT_H

#include <looper/measurement.h>
#include <alps/alea.h>
#include <cmath>
#include <complex>
#include <string>

struct gap_estimator : public looper::base_estimator
{
  template<typename T>
  static void initialize(T& m, bool is_bipartite, bool is_signed,
                         bool use_improved_estimator)
  {
    if (is_bipartite)
      looper::add_measurement(m, "Staggered Susceptibility [w=2pi/beta]",
                              is_signed);
    if (use_improved_estimator)
      looper::add_measurement(m, "Generalized Susceptibility [w=2pi/beta]",
                              is_signed);
  }

  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in)
  {
    if (m_in.has("Inverse Temperature")) {
      double beta = alps::RealObsevaluator(m_in["Inverse Temperature"]).mean();
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
  };

  struct estimate : public looper::base_estimator::estimate
  {
    std::complex<double> p;
    estimate() : p(0, 0) {}
    template<typename G, typename BIPARTITE>
    void start1(G const&, BIPARTITE, double t, int, int)
    { p -= looper::ctime(t); }
    template<typename G, typename BIPARTITE>
    void start1(G const&, BIPARTITE,
      looper::imaginary_time<boost::mpl::true_> const& t, int, int)
    { p -= t.ctime_; }
    template<typename G, typename BIPARTITE>
    void start2(G const& g, BIPARTITE, double t, int s, int c)
    { start1(g, BIPARTITE(), t, s, c); }
    template<typename G, typename BIPARTITE>
    void start2(G const& g, BIPARTITE,
      looper::imaginary_time<boost::mpl::true_> const& t, int s , int c)
    { start1(g, BIPARTITE(), t, s, c); }
    template<typename G, typename BIPARTITE>
    void term1(G const&, BIPARTITE, double t, int, int)
    { p += looper::ctime(t); }
    template<typename G, typename BIPARTITE>
    void term1(G const&, BIPARTITE,
      looper::imaginary_time<boost::mpl::true_> const& t, int, int)
    { p += t.ctime_; }
    template<typename G, typename BIPARTITE>
    void term2(G const& g, BIPARTITE, double t, int s, int c)
    { term1(g, BIPARTITE(), t, s, c); }
    template<typename G, typename BIPARTITE>
    void term2(G const& g, BIPARTITE,
      looper::imaginary_time<boost::mpl::true_> const& t, int s, int c)
    { term1(g, BIPARTITE(), t, s, c); }
  };

  template<typename QMC, typename BIPARTITE, typename IMPROVE>
  struct collector
  {
    double p;
    collector() : p(0) {}
    template<typename EST>
    collector operator+(EST const& cm)
    {
      p += looper::power2(cm.p);
      return *this;
    }
    template<typename M>
    void commit(M& m, double beta, int nrs, int /* nop */, double sign) const
    {
      using looper::power2;
      if (IMPROVE()) {
        if (BIPARTITE() && typename looper::is_path_integral<QMC>::type())
          m["Generalized Susceptibility [w=2pi/beta]"] <<
            sign * beta * p / power2(4*M_PI) / nrs;
      }
    }
  };

  template<typename QMC, typename BIPARTITE, typename IMPROVE>
  struct normal_estimator
  {
    template<typename M, typename G, typename OP>
    static void measure(M& m, G const& vg,
                        double beta, int nrs, int, double sign,
                        std::vector<int> const& spins,
                        std::vector<OP> const& operators,
                        std::vector<int>& spins_c)
    {
      if (!typename looper::is_path_integral<QMC>::type()) return;
      if (!BIPARTITE()) return;

      using looper::power2;
      typedef G graph_type;
      typedef typename alps::graph_traits<G>::site_iterator site_iterator;
      typedef typename std::vector<OP>::const_iterator operator_iterator;

      double smag = 0;
      site_iterator si, si_end;
      for (boost::tie(si, si_end) = sites(vg); si != si_end; ++si)
        smag += (0.5-spins[*si]) * looper::gauge(vg, *si);
      std::complex<double> smag_a(0, 0);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (operator_iterator oi = operators.begin(); oi != operators.end();
           ++oi) {
        if (oi->is_offdiagonal()) {
          std::complex<double> p = looper::ctime(oi->time());
          smag_a += p * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            smag += looper::gauge(vg, s) * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = source(bond(oi->pos(), vg), vg);
            unsigned int s1 = target(bond(oi->pos(), vg), vg);
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            smag += looper::gauge(vg, s0) * (1-2*spins_c[s0])
              + looper::gauge(vg, s1) * (1-2*spins_c[s1]);
          }
          smag_a -= p * smag;
        }
      }
      m["Staggered Susceptibility [w=2pi/beta]"] <<
        sign * beta * power2(smag_a) / power2(2*M_PI) / nrs;
    }
  };
};

#endif // GAP_MEASUREMENT_H
