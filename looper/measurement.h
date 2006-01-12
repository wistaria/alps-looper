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

#include "type.h"
#include "util.h"
#include <alps/alea.h>
#include <string>

namespace looper {

namespace measurement {

inline void add_measurement(alps::ObservableSet& m, std::string const& name,
                            bool is_signed = false)
{
  if (!m.has(name))
    m << make_observable(alps::RealObservable(name), is_signed);
}

inline void remove_measurement(alps::ObservableSet& m, std::string const& name)
{ if (m.has(name)) m.removeObservable(name); }

// for path integral
template<typename OP>
static inline void proceed(boost::mpl::true_, double& t, OP const& op)
{ t = op.time(); }
template<typename OP>
static inline void proceed(boost::mpl::false_, double&, OP const&) {}

// for sse
static inline void proceed(boost::mpl::true_, double& t) { t += 1; }
static inline void proceed(boost::mpl::false_, double&) {}

} // end namespace measurement


struct empty_improved_estimator
{
public:
  static void init(alps::ObservableSet&, bool, bool) {}
  static void evaluate(alps::ObservableSet&, double, int,
                       alps::ObservableSet const&) {}
  struct estimate
  {
    template<typename G, typename BIPARTITE>
    void at_zero(G const&, BIPARTITE, int, int) const {}
    template<typename G, typename BIPARTITE>
    void start(G const&, BIPARTITE, double, int, int) const {}
    template<typename G, typename BIPARTITE>
    void term(G const&, BIPARTITE, double, int, int) const {}
  };
  template<typename BIPARTITE>
  struct collector
  {
    template<typename EST>
    collector operator+(EST const&) const { return *this; }
    template<typename QMC>
    void commit(alps::ObservableSet&, QMC, double, int, int, double) const {}
  };
};

struct empty_normal_estimator
{
public:
  static void init(alps::ObservableSet&, bool, bool) {}
  static void evaluate(alps::ObservableSet&, double, int,
                       alps::ObservableSet const&) {}
  template<typename QMC, class G, class BIPARTITE, class OP>
  static void
  do_measurement(alps::ObservableSet const&, QMC, G const&, BIPARTITE,
                 double, int, int, double, std::vector<int> const&,
                 std::vector<OP> const&, std::vector<int> const&) {}
};


template<typename ESTIMATOR0, typename ESTIMATOR1>
struct improved_estimator_adaptor
{
public:
  typedef ESTIMATOR0 estimator0;
  typedef ESTIMATOR1 estimator1;
  static void init(alps::ObservableSet& m, bool is_bipartite, bool is_signed)
  {
    estimator0::init(m, is_bipartite, is_signed);
    estimator1::init(m, is_bipartite, is_signed);
  }
  static void evaluate(alps::ObservableSet& m_out, double beta, int nrs,
                       alps::ObservableSet const& m_in)
  {
    estimator0::evaluate(m_out, beta, nrs, m_in);
    estimator1::evaluate(m_out, beta, nrs, m_in);
  }

public:
  struct estimate : public estimator0::estimate, public estimator1::estimate
  {
    estimate() : estimator0::estimate(), estimator1::estimate() {}
    template<typename G, typename BIPARTITE>
    void at_zero(G const& g, BIPARTITE, int s, int c)
    {
      estimator0::estimate::at_zero(g, BIPARTITE(), s, c);
      estimator1::estimate::at_zero(g, BIPARTITE(), s, c);
    }
    template<typename G, typename BIPARTITE>
    void start(G const& g, BIPARTITE, double t, int s, int c)
    {
      estimator0::estimate::start(g, BIPARTITE(), t, s, c);
      estimator1::estimate::start(g, BIPARTITE(), t, s, c);
    }
    template<typename G, typename BIPARTITE>
    void term(G const& g, BIPARTITE, double t, int s, int c)
    {
      estimator0::estimate::term(g, BIPARTITE(), t, s, c);
      estimator1::estimate::term(g, BIPARTITE(), t, s, c);
    }
  };


private:
  template<typename EST, typename G, typename F, typename BIPARTITE,
           typename IMPROVE>
  struct accumulator_base
  {
    typedef EST estimate_t;
    typedef G lattice_graph_t;
    typedef F cluster_fragment_t;
    accumulator_base(std::vector<estimate_t> const&,
                     std::vector<cluster_fragment_t> const&,
                     lattice_graph_t const&) {}
    void at_zero(int, int, int) const {}
    void start(int, double, int, int) const {}
    void term(int, double, int, int) const {}
  };
  template<typename EST, typename G, typename F, typename BIPARTITE>
  struct accumulator_base<EST, G, F, BIPARTITE,
                          /* IMPROVE = */ boost::mpl::true_>
  {
    typedef EST estimate_t;
    typedef G lattice_graph_t;
    typedef F cluster_fragment_t;
    accumulator_base(std::vector<estimate_t>& es,
                    std::vector<cluster_fragment_t> const& fr,
                    lattice_graph_t const& vg)
      : estimates(es), fragments(fr), vgraph(vg) {}
    void at_zero(int p, int s, int c)
    { estimates[fragments[p].id].at_zero(vgraph, BIPARTITE(), s, c); }
    void start(int p, double t, int s, int c)
    { estimates[fragments[p].id].start(vgraph, BIPARTITE(), t, s, c); }
    void term(int p, double t, int s, int c)
    { estimates[fragments[p].id].term(vgraph, BIPARTITE(), t, s, c); }
    std::vector<estimate_t>& estimates;
    std::vector<cluster_fragment_t> const& fragments;
    lattice_graph_t const& vgraph;
  };

public:
  template<typename G, typename F, typename BIPARTITE, typename IMPROVE>
  struct accumulator :
    public accumulator_base<estimate, G, F, BIPARTITE, IMPROVE>
  {
    typedef accumulator_base<estimate, G, F, BIPARTITE, IMPROVE> super_type;
    accumulator(std::vector<estimate>& es, std::vector<F> const& fr,
                G const& vg)
      : accumulator_base<estimate, G, F, BIPARTITE, IMPROVE>(es, fr, vg) {}
  };

public:
  template<typename BIPARTITE>
  struct collector : public estimator0::template collector<BIPARTITE>,
                     public estimator1::template collector<BIPARTITE>
  {
    typedef typename estimator0::template collector<BIPARTITE> base0;
    typedef typename estimator1::template collector<BIPARTITE> base1;
    collector() : base0(), base1() {}
    template<typename EST>
    collector operator+(EST const& cm)
    {
      base0::operator+(cm);
      base1::operator+(cm);
      return *this;
    }
    template<typename QMC>
    void commit(alps::ObservableSet& m, QMC, double beta, int nrs, int nop,
                double sign) const
    {
      base0::commit(m, QMC(), beta, nrs, nop, sign);
      base1::commit(m, QMC(), beta, nrs, nop, sign);
    }
  };
};

template<typename ESTIMATOR0, typename ESTIMATOR1>
struct normal_estimator_adaptor
{
public:
  typedef ESTIMATOR0 estimator0;
  typedef ESTIMATOR1 estimator1;

  static void init(alps::ObservableSet& m, bool is_bipartite, bool is_signed)
  {
    estimator0::init(m, is_bipartite, is_signed);
    estimator1::init(m, is_bipartite, is_signed);
  }

  static void evaluate(alps::ObservableSet& m, double beta, int nrs,
                       alps::ObservableSet const& m_in)
  {
    estimator0::evaluate(m, beta, nrs, m_in);
    estimator1::evaluate(m, beta, nrs, m_in);
  }

  template<typename QMC, class G, class BIPARTITE, class OP>
  static void
  do_measurement(alps::ObservableSet& m, QMC, G const& vg, BIPARTITE,
                 double beta, int nrs, int nop, double sign,
                 std::vector<int> const& spins,
                 std::vector<OP> const& operators, std::vector<int>& spins_c)
  {
    estimator0::do_measurement(m, QMC(), vg, BIPARTITE(), beta, nrs, nop, sign,
                               spins, operators, spins_c);
    estimator1::do_measurement(m, QMC(), vg, BIPARTITE(), beta, nrs, nop, sign,
                               spins, operators, spins_c);
  }
};


struct susceptibility
{
  struct estimator_base
  {
    static void init(alps::ObservableSet& m, bool is_bipartite, bool is_signed,
                     bool use_improved_estimator)
    {
      using looper::measurement::add_measurement;
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

    static void evaluate(alps::ObservableSet& m,
                         double /* beta */, int /* nrs */,
                         alps::ObservableSet const& m_in)
    {
      using looper::measurement::remove_measurement;
      if (m_in.has("Magnetization^2") &&
          m_in.has("Magnetization^4")) {
        remove_measurement(m, "Binder Ratio of Magnetization");
        alps::RealObsevaluator obse_m2 = m_in["Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Magnetization^4"];
        alps::RealObsevaluator
          eval("Binder Ratio of Magnetization");
        eval = looper::power2(obse_m2) / obse_m4;
        m << eval;
      }
      if (m_in.has("Staggered Magnetization^2") &&
          m_in.has("Staggered Magnetization^4")) {
        remove_measurement(m,
          "Binder Ratio of Staggered Magnetization");
        alps::RealObsevaluator obse_m2 =
          m_in["Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 =
          m_in["Staggered Magnetization^4"];
        alps::RealObsevaluator
          eval("Binder Ratio of Staggered Magnetization");
        eval = looper::power2(obse_m2) / obse_m4;
        m << eval;
      }
      if (m_in.has("Generalized Magnetization^2") &&
          m_in.has("Generalized Magnetization^4")) {
        remove_measurement(m, "Binder Ratio of Generalized Magnetization");
        alps::RealObsevaluator obse_m2 = m_in["Generalized Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Generalized Magnetization^4"];
        alps::RealObsevaluator
          eval("Binder Ratio of Generalized Magnetization");
        eval = looper::power2(obse_m2) / obse_m4;
        m << eval;
      }
      if (m_in.has("Generalized Staggered Magnetization^2") &&
          m_in.has("Generalized Staggered Magnetization^4")) {
        remove_measurement(m,
          "Binder Ratio of Generalized Staggered Magnetization");
        alps::RealObsevaluator obse_m2 =
          m_in["Generalized Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 =
          m_in["Generalized Staggered Magnetization^4"];
        alps::RealObsevaluator
          eval("Binder Ratio of Generalized Staggered Magnetization");
        eval = looper::power2(obse_m2) / obse_m4;
        m << eval;
      }
    }
  };

  struct improved_estimator
  {
    static void init(alps::ObservableSet& m, bool is_bipartite, bool is_signed)
    { estimator_base::init(m, is_bipartite, is_signed, true); }

    static void evaluate(alps::ObservableSet& m, double beta, int nrs,
                         alps::ObservableSet const& m_in)
    { estimator_base::evaluate(m, beta, nrs, m_in); }

    struct estimate
    {
      double usize0, umag0, usize, umag;
      double ssize0, smag0, ssize, smag;
      estimate() : usize0(0), umag0(0), usize(0), umag(0),
                   ssize0(0), smag0(0), ssize(0), smag(0) {}
      template<typename G, typename BIPARTITE>
      void at_zero(G const& g, BIPARTITE, int s, int c)
      {
        usize0 += 0.5;
        umag0  += (0.5-c);
        if (BIPARTITE()) {
          double gg = gauge(g, s);
          ssize0 += gg * 0.5;
          smag0  += gg * (0.5-c);
        }
      }
      template<typename G, typename BIPARTITE>
      void start(G const& g, BIPARTITE, double t, int s, int c)
      {
        usize -= t * 0.5;
        umag  -= t * (0.5-c);
        if (BIPARTITE()) {
          double gg = gauge(g, s);
          ssize -= gg * t * 0.5;
          smag  -= gg * t * (0.5-c);
        }
      }
      template<typename G, typename BIPARTITE>
      void term(G const& g, BIPARTITE, double t, int s, int c)
      {
        usize += t * 0.5;
        umag  += t * (0.5-c);
        if (BIPARTITE()) {
          double gg = gauge(g, s);
          ssize += gg * t * 0.5;
          smag  += gg * t * (0.5-c);
        }
      }
    };

    template<typename EST, typename G, typename F, typename BIPARTITE,
             typename IMPROVE>
    struct accumulator_base
    {
      typedef EST estimate_t;
      typedef G lattice_graph_t;
      typedef F cluster_fragment_t;
      accumulator_base(std::vector<estimate_t> const&,
                       std::vector<cluster_fragment_t> const&,
                       lattice_graph_t const&) {}
      void at_zero(int, int, int) const {}
      void start(int, double, int, int) const {}
      void term(int, double, int, int) const {}
    };

    template<typename EST, typename G, typename F, typename BIPARTITE>
    struct accumulator_base<EST, G, F, BIPARTITE,
                            /* IMPROVE = */ boost::mpl::true_>
    {
      typedef EST estimate_t;
      typedef G lattice_graph_t;
      typedef F cluster_fragment_t;
      accumulator_base(std::vector<estimate_t>& es,
                       std::vector<cluster_fragment_t> const& fr,
                       lattice_graph_t const& vg)
        : estimates(es), fragments(fr), vgraph(vg) {}
      void at_zero(int p, int s, int c)
      { estimates[fragments[p].id].at_zero(vgraph, BIPARTITE(), s, c); }
      void start(int p, double t, int s, int c)
      { estimates[fragments[p].id].start(vgraph, BIPARTITE(), t, s, c); }
      void term(int p, double t, int s, int c)
      { estimates[fragments[p].id].term(vgraph, BIPARTITE(), t, s, c); }
      std::vector<estimate_t>& estimates;
      std::vector<cluster_fragment_t> const& fragments;
      lattice_graph_t const& vgraph;
    };

    template<typename G, typename F, typename BIPARTITE, typename IMPROVE>
    struct accumulator :
      public accumulator_base<estimate, G, F, BIPARTITE, IMPROVE>
    {
      typedef accumulator_base<estimate, G, F, BIPARTITE, IMPROVE> super_type;
      accumulator(std::vector<estimate>& es, std::vector<F> const& fr,
                  G const& vg)
        : accumulator_base<estimate, G, F, BIPARTITE, IMPROVE>(es, fr, vg) {}
    };

    template<typename BIPARTITE>
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
        if (BIPARTITE()) {
          ssize2 += power2(cm.ssize0);
          smag2  += power2(cm.smag0);
          ssize4 += power4(cm.ssize0);
          smag4  += power4(cm.smag0);
          ssize  += power2(cm.ssize);
          smag   += power2(cm.smag);
        }
        return *this;
      }
      template<typename QMC>
      void commit(alps::ObservableSet& m, QMC, double beta, int nrs, int nop,
                  double sign) const
      {
        m["Magnetization"] << 0.0;
        m["Magnetization Density"] << 0.0;
        m["Magnetization^2"] << sign * umag2;
        m["Magnetization^4"] << sign * (3 * umag2 * umag2 - 2 * umag4);
        m["Susceptibility"]
          << (typename is_path_integral<QMC>::type() ?
              sign * beta * umag / nrs :
              sign * beta * (dip(umag, nop) + umag2) / (nop + 1) / nrs);
        m["Generalized Magnetization^2"] << sign * usize2;
        m["Generalized Magnetization^4"]
          << sign * (3 * usize2 * usize2 - 2 * usize4);
        m["Generalized Susceptibility"]
          << (typename is_path_integral<QMC>::type() ?
              sign * beta * usize / nrs :
              sign * beta * (dip(usize, nop) + usize2) / (nop + 1) / nrs);
        if (BIPARTITE()) {
          m["Staggered Magnetization"] << 0.0;
          m["Staggered Magnetization Density"] << 0.0;
          m["Staggered Magnetization^2"] << sign * smag2;
          m["Staggered Magnetization^4"]
            << sign * (3 * smag2 * smag2 - 2 * smag4);
          m["Staggered Susceptibility"]
            << (typename is_path_integral<QMC>::type() ?
                sign * beta * smag /nrs :
                sign * beta * (dip(smag, nop) + smag2) / (nop + 1) / nrs);
          m["Generalized Staggered Magnetization^2"] << sign * ssize2;
          m["Generalized Staggered Magnetization^4"]
            << sign * (3 * ssize2 * ssize2 - 2 * ssize4);
          m["Generalized Staggered Susceptibility"]
            << (typename is_path_integral<QMC>::type() ?
                sign * beta * ssize / nrs :
                sign * beta * (dip(ssize, nop) + ssize2) / (nop + 1) / nrs);
        }
      }
    };
  };

  struct normal_estimator
  {
    static void init(alps::ObservableSet& m, bool is_bipartite, bool is_signed)
    { estimator_base::init(m, is_bipartite, is_signed, false); }

    static void evaluate(alps::ObservableSet& m, double beta, int nrs,
                         alps::ObservableSet const& m_in)
    { estimator_base::evaluate(m, beta, nrs, m_in); }

    template<typename QMC, class G, class BIPARTITE, class OP>
    static void
    do_measurement(alps::ObservableSet& m, QMC, G const& vg, BIPARTITE,
                   double beta, int nrs, int nop, double sign,
                   std::vector<int> const& spins,
                   std::vector<OP> const& operators, std::vector<int>& spins_c)
    {
      typedef G graph_type;
      typedef typename alps::graph_traits<G>::site_iterator site_iterator;
      typedef typename std::vector<OP>::const_iterator operator_iterator;

      double umag = 0;
      double smag = 0;
      site_iterator si, si_end;
      for (boost::tie(si, si_end) = sites(vg); si != si_end; ++si) {
        umag += 0.5-spins[*si];
        if (BIPARTITE()) smag += (0.5-spins[*si]) * gauge(vg, *si);
      }
      m["Magnetization"] << sign * umag;
      m["Magnetization Density"] << sign * umag / nrs;
      m["Magnetization^2"] << sign * power2(umag);
      m["Magnetization^4"] << sign * power4(umag);
      if (BIPARTITE()) {
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
          if (BIPARTITE()) smag_a += t * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += 1-2*spins_c[s];
            if (BIPARTITE()) smag += gauge(vg, s) * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = source(bond(oi->pos(), vg), vg);
            unsigned int s1 = target(bond(oi->pos(), vg), vg);
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            umag += 1-2*spins_c[s0] + 1-2*spins_c[s1];
            if (BIPARTITE()) smag += gauge(vg, s0) * (1-2*spins_c[s0])
                               + gauge(vg, s1) * (1-2*spins_c[s1]);
          }
          umag_a -= t * umag;
          if (BIPARTITE()) smag_a -= t * smag;
        }
        measurement::proceed(typename is_sse<QMC>::type(), t);
      }
      if (typename is_path_integral<QMC>::type()) {
        umag_a += beta * umag;
        m["Susceptibility"] << sign * beta * power2(umag_a) / nrs;
        if (BIPARTITE()) {
          smag_a += beta * smag;
          m["Staggered Susceptibility"] << sign * beta * power2(smag_a) / nrs;
        }
      } else {
        umag_a += nop * umag;
        m["Susceptibility"]
          << sign * beta * (dip(power2(umag_a), nop) + power2(umag))
          / (nop + 1) / nrs;
        if (BIPARTITE()) {
          smag_a += nop * smag;
          m["Staggered Susceptibility"]
            << sign * beta * (dip(power2(smag_a), nop) + power2(smag))
            / (nop + 1) / nrs;
        }
      }
    }
  };
};

typedef susceptibility::improved_estimator improved_estimator;
typedef susceptibility::normal_estimator   normal_estimator;

} // end namespace looper

#endif // LOOPER_MEASUREMENT_H
