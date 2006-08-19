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

template<typename VLAT, typename TIME>
struct gap_estimator
{
  typedef VLAT virtual_lattice_t;
  typedef TIME time_t;

  template<typename M>
  void initialize(M& m, alps::Parameters const& /* params */,
                  virtual_lattice_t const& /* vlat */,
                  bool is_bipartite, bool is_signed,
                  bool use_improved_estimator)
  {
    if (is_bipartite)
      looper::add_measurement(m, "Staggered Susceptibility [w=2pi/beta]",
                              is_signed);
    if (use_improved_estimator)
      looper::add_measurement(m, "Generalized Susceptibility [w=2pi/beta]",
                              is_signed);
  }

  static void evaluate(alps::ObservableSet& m,
                       alps::Parameters const& /* params */,
                       alps::ObservableSet const& m_in)
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

  // improved estimator

  struct estimate
  {
    std::complex<double> p;
    void init() { p = std::complex<double>(0,0); }
    void start_s(virtual_lattice_t const&, double t, int, int)
    { p -= looper::ctime(t); }
    void start_s(virtual_lattice_t const&,
      looper::imaginary_time<boost::mpl::true_> const& t, int, int)
    { p -= t.ctime_; }
    void start_b(virtual_lattice_t const& vlat, double t, int, int s, int c)
    { start_s(vlat, t, s, c); }
    void start_b(virtual_lattice_t const& vlat,
      looper::imaginary_time<boost::mpl::true_> const& t, int, int s , int c)
    { start_s(vlat, t, s, c); }
    void term_s(virtual_lattice_t const&, double t, int, int)
    { p += looper::ctime(t); }
    void term_s(virtual_lattice_t const&,
      looper::imaginary_time<boost::mpl::true_> const& t, int, int)
    { p += t.ctime_; }
    void term_b(virtual_lattice_t const& vlat, double t, int, int s, int c)
    { term_s(vlat, t, s, c); }
    void term_b(virtual_lattice_t const& vlat,
      looper::imaginary_time<boost::mpl::true_> const& t, int, int s, int c)
    { term_s(vlat, t, s, c); }
    void at_bot(virtual_lattice_t const& vlat, double t, int s, int c)
    { start_s(vlat, t, s, c); }
    void at_bot(virtual_lattice_t const& vlat,
      looper::imaginary_time<boost::mpl::true_> const& t, int s, int c)
    { start_s(vlat, t, s, c); }
    void at_top(virtual_lattice_t const& vlat, double t, int s, int c)
    { term_s(vlat, t, s, c); }
    void at_top(virtual_lattice_t const& vlat,
      looper::imaginary_time<boost::mpl::true_> const& t, int s, int c)
    { term_s(vlat, t, s, c); }
  };
  void init_estimate(estimate& est) const { est.init(); }

  template<typename QMC, typename IMPROVE>
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
    void commit(M& m, virtual_lattice_t const& vlat, bool is_bipartite,
                double beta, int, double sign) const
    {
      if (is_bipartite)
        m["Generalized Susceptibility [w=2pi/beta]"] <<
          sign * beta * p / looper::power2(4*M_PI) / num_sites(vlat.rgraph());
    }
  };

  // normal estimator

  template<typename QMC, typename M, typename OP>
  void measure(M& m, virtual_lattice_t const vlat,
               bool is_bipartite, bool /* use_improved_estimator */,
               double beta, double sign,
               std::vector<int> const& spins,
               std::vector<OP> const& operators,
               std::vector<int>& spins_c)
  {
    if (!typename looper::is_path_integral<QMC>::type()) return;
    if (!is_bipartite) return;

    int nrs = num_sites(vlat.rgraph());

    double smag = 0;
    typename virtual_lattice_t::virtual_site_iterator
      si, si_end;
    for (boost::tie(si, si_end) = vsites(vlat); si != si_end; ++si)
      smag += (0.5-spins[*si]) * looper::gauge(vlat, *si);
    std::complex<double> smag_a(0, 0);
    std::copy(spins.begin(), spins.end(), spins_c.begin());
    for (typename std::vector<OP>::const_iterator oi = operators.begin();
         oi != operators.end(); ++oi) {
      if (oi->is_offdiagonal()) {
        std::complex<double> p = looper::ctime(oi->time());
        smag_a += p * smag;
        if (oi->is_site()) {
          unsigned int s = oi->pos();
          spins_c[s] ^= 1;
          smag += looper::gauge(vlat, s) * (1-2*spins_c[s]);
        } else {
          unsigned int s0 = vsource(oi->pos(), vlat);
          unsigned int s1 = vtarget(oi->pos(), vlat);
          spins_c[s0] ^= 1;
          spins_c[s1] ^= 1;
          smag += looper::gauge(vlat, s0) * (1-2*spins_c[s0])
            + looper::gauge(vlat, s1) * (1-2*spins_c[s1]);
        }
        smag_a -= p * smag;
      }
    }
    m["Staggered Susceptibility [w=2pi/beta]"] <<
      sign * beta * looper::power2(smag_a) / looper::power2(2*M_PI) / nrs;
  }
};

#endif // GAP_MEASUREMENT_H
