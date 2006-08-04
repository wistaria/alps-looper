/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_MEASUREMENT_H
#define LOOPER_MEASUREMENT_H

#include "divide_if_positive.h"
#include "lattice.h"
#include "power.h"
#include "type.h"

#include <alps/fixed_capacity_vector.h>
#include <alps/lattice/graph_traits.h>
#include <alps/alea.h>
#include <boost/call_traits.hpp>
#include <boost/mpl/bool.hpp>
#include <string>
#include <vector>

namespace looper {

//
// helper functions
//

inline void add_measurement(alps::ObservableSet& m, std::string const& name,
                            bool is_signed = false)
{
  if (!m.has(name))
    m << make_observable(alps::RealObservable(name), is_signed);
}

namespace measurement {

// for path integral
template<typename OP>
inline void proceed(boost::mpl::true_, double& t, OP const& op)
{ t = op.time(); }
template<typename OP>
inline void proceed(boost::mpl::false_, double&, OP const&) {}

// for sse
inline void proceed(boost::mpl::true_, double& t) { t += 1; }
inline void proceed(boost::mpl::false_, double&) {}


//
// traits classes
//

template<typename ESTIMATOR>
struct estimate
{
  typedef typename ESTIMATOR::estimate type;
};

template<typename ESTIMATOR, typename G, typename T, typename F,
  typename IMPROVE>
struct accumulator
{
  typedef typename estimate<ESTIMATOR>::type estimate_t;
  typedef typename
    ESTIMATOR::template accumulator<estimate_t, G, T, F, IMPROVE>
    type;
};

template<typename ESTIMATOR, typename QMC, typename IMPROVE>
struct collector
{
  typedef typename ESTIMATOR::template collector<QMC, IMPROVE> type;
};

template<typename ESTIMATOR, typename QMC, typename IMPROVE>
struct normal_estimator
{
  typedef typename
    ESTIMATOR::template normal_estimator<QMC, IMPROVE>
    type;
};

} // end namespace measurement


//
// base_estimator
//

struct base_estimator
{
  template<typename M>
  static void initialize(M& /* m */,
                         bool /* is_bipartite */,
                         bool /* is_signed */,
                         bool /* use_improved_estimator */) {}
  static void evaluate(alps::ObservableSet&, alps::ObservableSet const&) {}

  // improved estimator

  struct estimate
  {
    template<typename G>
    void init(G const&) {}
    template<typename G>
    void start_s(G const&, double, int, int) const {}
    template<typename G>
    void start_b(G const&, double, int, int, int) const {}
    template<typename G>
    void term_s(G const&, double, int, int) const {}
    template<typename G>
    void term_b(G const&, double, int, int, int) const {}
    template<typename G>
    void at_bot(G const&, double, int, int) const {}
    template<typename G>
    void at_top(G const&, double, int, int) const {}
  };
  template<typename EST, typename G, typename T, typename F,
           typename IMPROVE>
  struct accumulator
  {
    typedef EST estimate_t;
    typedef G lattice_graph_t;
    typedef T time_t;
    typedef typename boost::call_traits<time_t>::param_type time_pt;
    typedef F cluster_fragment_t;
    accumulator(int, std::vector<estimate_t> const&,
                std::vector<cluster_fragment_t> const&,
                virtual_lattice<lattice_graph_t> const&) {}
    void start_s(int, time_pt, int, int) const {}
    void start_b(int, int, time_pt, int, int, int, int, int) const {}
    void term_s(int, time_pt, int, int) const {}
    void term_b(int, int, time_pt, int, int, int, int, int) const {}
    void at_bot(int, time_pt, int, int) const {}
    void at_top(int, time_pt, int, int) const {}
  };
  template<typename EST, typename G, typename T, typename F>
  struct accumulator<EST, G, T, F, boost::mpl::true_>
  {
    typedef EST estimate_t;
    typedef G lattice_graph_t;
    typedef T time_t;
    typedef typename boost::call_traits<time_t>::param_type time_pt;
    typedef F cluster_fragment_t;
    accumulator(int nc, std::vector<estimate_t>& es,
                std::vector<cluster_fragment_t> const& fr,
                virtual_lattice<lattice_graph_t> const& vl)
      : estimates(es), fragments(fr), vgraph(vl.graph())
    {
      estimate_t e;
      e.init(vl);
      estimates.resize(0); estimates.resize(nc, e);
    }
    void start_s(int p, time_pt t, int s, int c)
    { estimates[fragments[p].id].start_s(vgraph, t, s, c); }
    void start_b(int p0, int p1, time_pt t, int b, int s0, int s1,
                 int c0, int c1)
    {
      estimates[fragments[p0].id].start_b(vgraph, t, b, s0, c0);
      estimates[fragments[p1].id].start_b(vgraph, t, b, s1, c1);
    }
    void term_s(int p, time_pt t, int s, int c)
    { estimates[fragments[p].id].term_s(vgraph, t, s, c); }
    void term_b(int p0, int p1, time_pt t, int b, int s0, int s1,
                int c0, int c1)
    {
      estimates[fragments[p0].id].term_b(vgraph, t, b, s0, c0);
      estimates[fragments[p1].id].term_b(vgraph, t, b, s1, c1);
    }
    void at_bot(int p, time_pt t, int s, int c)
    { estimates[fragments[p].id].at_bot(vgraph, t, s, c); }
    void at_top(int p, time_pt t, int s, int c)
    { estimates[fragments[p].id].at_top(vgraph, t, s, c); }
    std::vector<estimate_t>& estimates;
    std::vector<cluster_fragment_t> const& fragments;
    lattice_graph_t const& vgraph;
  };
  template<typename QMC, typename IMPROVE>
  struct collector
  {
    template<typename EST>
    collector operator+(EST const&) const { return *this; }
    template<typename M>
    void commit(M&, bool, double, int, int, double) const {}
  };

  // normal estimator

  template<typename QMC, typename IMPROVE>
  struct normal_estimator
  {
    template<typename M, typename G, typename OP>
    static void measure(M const&, bool, G const&, virtual_lattice<G> const&,
                        double, int, int, double,
                        std::vector<int> const&, std::vector<OP> const&,
                        std::vector<int> const&) {}
  };
};


//
// estimator_adaptor
//

template<typename ESTIMATOR1, typename ESTIMATOR2>
struct composite_estimator : public base_estimator
{
  typedef ESTIMATOR1 estimator1;
  typedef ESTIMATOR2 estimator2;

  template<typename T>
  static void initialize(T& m, bool is_bipartite, bool is_signed,
                         bool use_improved_estimator)
  {
    estimator1::initialize(m, is_bipartite, is_signed, use_improved_estimator);
    estimator2::initialize(m, is_bipartite, is_signed, use_improved_estimator);
  };

  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in)
  {
    estimator1::evaluate(m, m_in);
    estimator2::evaluate(m, m_in);
  }

  struct estimate : public estimator1::estimate, public estimator2::estimate
  {
    template<typename G>
    void init(G const& g)
    {
      estimator1::estimate::init(g);
      estimator2::estimate::init(g);
    }
    template<typename G>
    void start_s(G const& g, double t, int s, int c)
    {
      estimator1::estimate::start_s(g, t, s, c);
      estimator2::estimate::start_s(g, t, s, c);
    }
    template<typename G, typename T>
    void start_s(G const& g, T const& t, int s, int c)
    {
      estimator1::estimate::start_s(g, t, s, c);
      estimator2::estimate::start_s(g, t, s, c);
    }
    template<typename G>
    void start_b(G const& g, double t, int b, int s, int c)
    {
      estimator1::estimate::start_b(g, t, b, s, c);
      estimator2::estimate::start_b(g, t, b, s, c);
    }
    template<typename G, typename T>
    void start_b(G const& g, T const& t, int b, int s, int c)
    {
      estimator1::estimate::start_b(g, t, b, s, c);
      estimator2::estimate::start_b(g, t, b, s, c);
    }
    template<typename G>
    void term_s(G const& g, double t, int s, int c)
    {
      estimator1::estimate::term_s(g, t, s, c);
      estimator2::estimate::term_s(g, t, s, c);
    }
    template<typename G, typename T>
    void term_s(G const& g, T const& t, int s, int c)
    {
      estimator1::estimate::term_s(g, t, s, c);
      estimator2::estimate::term_s(g, t, s, c);
    }
    template<typename G>
    void term_b(G const& g, double t, int b, int s, int c)
    {
      estimator1::estimate::term_b(g, t, b, s, c);
      estimator2::estimate::term_b(g, t, b, s, c);
    }
    template<typename G, typename T>
    void term_b(G const& g, T const& t, int b, int s, int c)
    {
      estimator1::estimate::term_b(g, t, b, s, c);
      estimator2::estimate::term_b(g, t, b, s, c);
    }
    template<typename G>
    void at_bot(G const& g, double t, int s, int c)
    {
      estimator1::estimate::at_bot(g, t, s, c);
      estimator2::estimate::at_bot(g, t, s, c);
    }
    template<typename G, typename T>
    void at_bot(G const& g, T const& t, int s, int c)
    {
      estimator1::estimate::at_bot(g, t, s, c);
      estimator2::estimate::at_bot(g, t, s, c);
    }
    template<typename G>
    void at_top(G const& g, double t, int s, int c)
    {
      estimator1::estimate::at_top(g, t, s, c);
      estimator2::estimate::at_top(g, t, s, c);
    }
    template<typename G, typename T>
    void at_top(G const& g, T const& t, int s, int c)
    {
      estimator1::estimate::at_top(g, t, s, c);
      estimator2::estimate::at_top(g, t, s, c);
    }
  };

  template<typename QMC, typename IMPROVE>
  struct collector
    : public estimator1::template collector<QMC, IMPROVE>,
      public estimator2::template collector<QMC, IMPROVE>
  {
    typedef typename estimator1::template collector<QMC, IMPROVE>
      base1;
    typedef typename estimator2::template collector<QMC, IMPROVE>
      base2;
    collector() : base1(), base2() {}
    template<typename EST>
    collector operator+(EST const& cm)
    {
      base1::operator+(cm);
      base2::operator+(cm);
      return *this;
    }
    template<typename M>
    void commit(M& m, bool is_bipartite, double beta, int nrs, int nop,
                double sign) const
    {
      base1::commit(m, is_bipartite, beta, nrs, nop, sign);
      base2::commit(m, is_bipartite, beta, nrs, nop, sign);
    }
  };

  template<typename QMC, typename IMPROVE>
  struct normal_estimator
  {
    template<typename M, typename G, typename OP>
    static void measure(M& m, bool is_bipartite,
                        G const& rg, virtual_lattice<G> const& vl,
                        double beta, int nrs, int nop, double sign,
                        std::vector<int> const& spins,
                        std::vector<OP> const& operators,
                        std::vector<int>& spins_c)
    {
      estimator1::template normal_estimator<QMC, IMPROVE>::
        measure(m, is_bipartite, rg, vl, beta, nrs, nop, sign, spins,
                operators, spins_c);
      estimator2::template normal_estimator<QMC, IMPROVE>::
        measure(m, is_bipartite, rg, vl, beta, nrs, nop, sign, spins,
                operators, spins_c);
    }
  };
};


//
// energy_estimator
//

struct energy_estimator
{
  template<typename T>
  static void initialize(T& m, bool is_signed)
  {
    add_measurement(m, "Energy", is_signed);
    add_measurement(m, "Energy Density", is_signed);
    add_measurement(m, "Energy^2", is_signed);
  }

  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in)
  {
    if (m_in.has("Inverse Temperature") && m_in.has("Number of Sites") &&
        m_in.has("Energy") && m_in.has("Energy^2")) {
      double beta = alps::RealObsevaluator(m_in["Inverse Temperature"]).mean();
      double nrs = alps::RealObsevaluator(m_in["Number of Sites"]).mean();
      alps::RealObsevaluator obse_e = m_in["Energy"];
      alps::RealObsevaluator obse_e2 = m_in["Energy^2"];
      alps::RealObsevaluator eval("Specific Heat");
      eval = power2(beta) * (obse_e2 - power2(obse_e)) / nrs;
      m.addObservable(eval);
    }
  }

  // no improved estimator

  // normal_estimator

  static void measure(alps::ObservableSet& m, double beta, int nrs, int nop,
                      double sign, double ene)
  {
    m["Energy"] << sign * ene;
    m["Energy Density"] << sign * ene / nrs;
    m["Energy^2"] << sign * (power2(ene) - nop / power2(beta));
  }
};

//
// susceptibility_estimator
//

struct susceptibility_estimator : public base_estimator
{
  template<typename T>
  static void initialize(T& m, bool is_bipartite, bool is_signed,
                         bool use_improved_estimator)
  {
    add_measurement(m, "Magnetization", is_signed);
    add_measurement(m, "Magnetization Density", is_signed);
    add_measurement(m, "Magnetization^2", is_signed);
    add_measurement(m, "Magnetization^4", is_signed);
    add_measurement(m, "Susceptibility", is_signed);
    if (use_improved_estimator) {
      add_measurement(m, "Generalized Magnetization^2", is_signed);
      add_measurement(m, "Generalized Magnetization^4", is_signed);
      add_measurement(m, "Generalized Susceptibility", is_signed);
    }
    if (is_bipartite) {
      add_measurement(m, "Staggered Magnetization", is_signed);
      add_measurement(m, "Staggered Magnetization Density", is_signed);
      add_measurement(m, "Staggered Magnetization^2", is_signed);
      add_measurement(m, "Staggered Magnetization^4", is_signed);
      add_measurement(m, "Staggered Susceptibility", is_signed);
      if (use_improved_estimator) {
        add_measurement(m, "Generalized Staggered Magnetization^2",
                        is_signed);
        add_measurement(m, "Generalized Staggered Magnetization^4",
                        is_signed);
        add_measurement(m, "Generalized Staggered Susceptibility",
                        is_signed);
      }
    }
  }

  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in)
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
      alps::RealObsevaluator eval("Binder Ratio of Generalized Magnetization");
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

  struct estimate
  {
    double usize0, umag0, usize, umag;
    double ssize0, smag0, ssize, smag;
    template<typename G>
    void init(G const&)
    {
      usize0 = 0;
      umag0 = 0;
      usize = 0;
      umag = 0;
      ssize0 = 0;
      smag0 = 0;
      ssize = 0;
      smag = 0;
    }
    template<typename G>
    void start_s(G const& g, double t, int s, int c)
    {
      usize -= t * 0.5;
      umag  -= t * (0.5-c);
      double gg = gauge(g, s);
      ssize -= gg * t * 0.5;
      smag  -= gg * t * (0.5-c);
    }
    template<typename G>
    void start_b(G const& g, double t, int, int s, int c)
    { start_s(g, t, s, c); }
    template<typename G>
    void term_s(G const& g, double t, int s, int c)
    {
      usize += t * 0.5;
      umag  += t * (0.5-c);
      double gg = gauge(g, s);
      ssize += gg * t * 0.5;
      smag  += gg * t * (0.5-c);
    }
    template<typename G>
    void term_b(G const& g, double t, int, int s, int c)
    { term_s(g, t, s, c); }
    template<typename G>
    void at_bot(G const& g, double t, int s, int c)
    {
      start_s(g, t, s, c);
      usize0 += 0.5;
      umag0  += (0.5-c);
      double gg = gauge(g, s);
      ssize0 += gg * 0.5;
      smag0  += gg * (0.5-c);
    }
    template<typename G>
    void at_top(G const& g, double t, int s, int c)
    { term_s(g, t, s, c); }
  };

  template<typename QMC, typename IMPROVE>
  struct collector
  {
    double usize2, umag2, usize4, umag4, usize, umag;
    double ssize2, smag2, ssize4, smag4, ssize, smag;
    collector() : usize2(0), umag2(0), usize4(0), umag4(0),
                  usize(0), umag(0),
                  ssize2(0), smag2(0), ssize4(0), smag4(0),
                  ssize(0), smag(0) {}
    template<typename EST>
    collector operator+(EST const& cm)
    {
      usize2 += power2(cm.usize0);
      umag2  += power2(cm.umag0);
      usize4 += power4(cm.usize0);
      umag4  += power4(cm.umag0);
      usize  += power2(cm.usize);
      umag   += power2(cm.umag);
      ssize2 += power2(cm.ssize0);
      smag2  += power2(cm.smag0);
      ssize4 += power4(cm.ssize0);
      smag4  += power4(cm.smag0);
      ssize  += power2(cm.ssize);
      smag   += power2(cm.smag);
      return *this;
    }
    template<typename M>
    void commit(M& m, bool is_bipartite, double beta, int nrs, int nop,
                double sign) const
    {
      m["Magnetization"] << 0.0;
      m["Magnetization Density"] << 0.0;
      m["Magnetization^2"] << sign * umag2;
      m["Magnetization^4"] << sign * (3 * umag2 * umag2 - 2 * umag4);
      m["Susceptibility"]
        << (typename is_sse<QMC>::type() ?
            sign * beta * (dip(umag, nop) + umag2) / (nop + 1) / nrs :
            sign * beta * umag / nrs);
      m["Generalized Magnetization^2"] << sign * usize2;
      m["Generalized Magnetization^4"]
        << sign * (3 * usize2 * usize2 - 2 * usize4);
      m["Generalized Susceptibility"]
        << (typename is_sse<QMC>::type() ?
            sign * beta * (dip(usize, nop) + usize2) / (nop + 1) / nrs :
            sign * beta * usize / nrs);
      if (is_bipartite) {
        m["Staggered Magnetization"] << 0.0;
        m["Staggered Magnetization Density"] << 0.0;
        m["Staggered Magnetization^2"] << sign * smag2;
        m["Staggered Magnetization^4"]
          << sign * (3 * smag2 * smag2 - 2 * smag4);
        m["Staggered Susceptibility"]
          << (typename is_sse<QMC>::type() ?
              sign * beta * (dip(smag, nop) + smag2) / (nop + 1) / nrs :
              sign * beta * smag /nrs);
        m["Generalized Staggered Magnetization^2"] << sign * ssize2;
        m["Generalized Staggered Magnetization^4"]
          << sign * (3 * ssize2 * ssize2 - 2 * ssize4);
        m["Generalized Staggered Susceptibility"]
          << (typename is_sse<QMC>::type() ?
              sign * beta * (dip(ssize, nop) + ssize2) / (nop + 1) / nrs :
              sign * beta * ssize / nrs);
      }
    }
  };

  template<typename QMC, typename IMPROVE>
  struct normal_estimator
  {
    template<typename M, typename G, typename OP>
    static void measure(M& m, bool is_bipartite,
                        G const& /* rg */, virtual_lattice<G> const& vl,
                        double beta, int nrs, int nop, double sign,
                        std::vector<int> const& spins,
                        std::vector<OP> const& operators,
                        std::vector<int>& spins_c)
    {
      if (IMPROVE()) return;

      typedef G graph_type;
      typedef typename alps::graph_traits<G>::site_iterator site_iterator;
      typedef typename std::vector<OP>::const_iterator operator_iterator;

      double umag = 0;
      double smag = 0;
      site_iterator si, si_end;
      for (boost::tie(si, si_end) = sites(vl); si != si_end; ++si) {
        umag += 0.5-spins[*si];
        smag += (0.5-spins[*si]) * gauge(vl, *si);
      }
      m["Magnetization"] << sign * umag;
      m["Magnetization Density"] << sign * umag / nrs;
      m["Magnetization^2"] << sign * power2(umag);
      m["Magnetization^4"] << sign * power4(umag);
      if (is_bipartite) {
        m["Staggered Magnetization"] << sign * smag;
        m["Staggered Magnetization Density"] << sign * smag / nrs;
        m["Staggered Magnetization^2"] << sign * power2(smag);
        m["Staggered Magnetization^4"] << sign * power4(smag);
      }
      double umag_a = 0; /* 0 * umag; */
      double smag_a = 0; /* 0 * smag; */
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      double t = 0;
      for (operator_iterator oi = operators.begin(); oi != operators.end();
           ++oi) {
        if (oi->is_offdiagonal()) {
          measurement::proceed(typename is_path_integral<QMC>::type(), t, *oi);
          umag_a += t * umag;
          smag_a += t * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += 1-2*spins_c[s];
            smag += gauge(vl, s) * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = vsource(oi->pos(), vl);
            unsigned int s1 = vtarget(oi->pos(), vl);
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            umag += 1-2*spins_c[s0] + 1-2*spins_c[s1];
            smag += gauge(vl, s0) * (1-2*spins_c[s0])
              + gauge(vl, s1) * (1-2*spins_c[s1]);
          }
          umag_a -= t * umag;
          smag_a -= t * smag;
        }
        measurement::proceed(typename is_sse<QMC>::type(), t);
      }
      if (typename is_path_integral<QMC>::type()) {
        umag_a += umag;
        m["Susceptibility"] << sign * beta * power2(umag_a) / nrs;
        smag_a += smag;
        if (is_bipartite)
          m["Staggered Susceptibility"] << sign * beta * power2(smag_a) / nrs;
      } else {
        umag_a += nop * umag;
        m["Susceptibility"]
          << sign * beta * (dip(power2(umag_a), nop) + power2(umag))
          / (nop + 1) / nrs;
        smag_a += nop * smag;
        if (is_bipartite)
          m["Staggered Susceptibility"]
            << sign * beta * (dip(power2(smag_a), nop) + power2(smag))
            / (nop + 1) / nrs;
      }
    }
  };
};


//
// stiffness_estimator
//

template<unsigned int MAX_DIM>
struct stiffness_estimator : public base_estimator
{
  template<typename M>
  static void initialize(M& m,
                         bool /* is_bipartite */,
                         bool is_signed,
                         bool /* use_improved_estimator */)
  {
    add_measurement(m, "Stiffness", is_signed);
  }
  // static void evaluate(alps::ObservableSet&, alps::ObservableSet const&) {}

  // improved estimator

  struct estimate
  {
    alps::fixed_capacity_vector<double, MAX_DIM> winding;
    template<typename G>
    void init(G const& vl)
    {
      int dim = boost::get_property(vl.graph(), dimension_t());
      if (dim > MAX_DIM) {
        std::cerr << "Spatial dimension (=" << dim << ") is too large.  "
                  << "Stiffness will be measured only for the first "
                  << MAX_DIM << " dimensions\n";
        dim = MAX_DIM;
      }
      winding.resize(dim);
    }
    template<typename G>
    void start_s(G const&, double, int, int) const {}
    template<typename G>
    void start_b(G const&, double, int, int, int) const {}
    template<typename G>
    void term_s(G const&, double, int, int) const {}
    template<typename G>
    void term_b(G const&, double, int, int, int) const {}
    template<typename G>
    void at_bot(G const&, double, int, int) const {}
    template<typename G>
    void at_top(G const&, double, int, int) const {}
  };

  template<typename QMC, typename IMPROVE>
  struct normal_estimator
  {
    template<typename M, typename G, typename OP>
    static void measure(M& m, bool is_bipartite,
                        G const& rg, virtual_lattice<G> const& vl,
                        double beta, int nrs, int nop, double sign,
                        std::vector<int> const& spins,
                        std::vector<OP> const& operators,
                        std::vector<int>& spins_c)
    {
      // if (IMPROVE()) return;

      typedef typename std::vector<OP>::const_iterator operator_iterator;

      int dim = boost::get_property(vl.graph(), dimension_t());
      std::valarray<double> winding(0., dim);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (operator_iterator oi = operators.begin(); oi != operators.end();
           ++oi) {
        if (oi->is_offdiagonal()) {
          if (oi->is_bond()) {
            double s = 1-2*spins_c[vsource(oi->pos(), vl)];
            alps::coordinate_type vr =
              get(bond_vector_relative_t(), vl.graph(), oi->pos());
            for (int i = 0; i < dim; ++i) winding[i] += s * vr[i];
            spins_c[vsource(oi->pos(), vl)] ^= 1;
            spins_c[vtarget(oi->pos(), vl)] ^= 1;
          } else {
            spins_c[oi->pos()] ^= 1;
          }
        }
      }

      double w = 0;
      for (int i = 0; i < dim; ++i) w += power2(winding[i]);
      m["Stiffness"] << sign * w / (beta * dim);
    }
  };
};

} // end namespace looper

#endif // LOOPER_MEASUREMENT_H
