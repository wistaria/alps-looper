/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOP_WORKER_H
#define LOOP_WORKER_H

#include <alps/alea.h>
#include <alps/scheduler.h>
#include <looper.h>

template<class MCRUN = alps::scheduler::LatticeModelMCRun<looper::graph_type> >
class qmc_worker_base : public MCRUN
{
public:
  typedef MCRUN                           super_type;
  typedef typename super_type::graph_type graph_type;
  typedef typename looper::graph_traits<graph_type>::site_iterator
    site_iterator;
  typedef typename looper::graph_traits<graph_type>::bond_iterator
    bond_iterator;

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n),
      mcs_therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
      mcs_sweep_(p["SWEEPS"]),
      vlat_(), gtab_(),
      r_graph(*super_type::engine_ptr, looper::random_choice<>()),
      mcs_(0)
  {
    //
    // setup model
    //

    looper::model_parameter mp(p, *this);
    bool is_signed = mp.is_signed();
    bool is_classically_frustrated = mp.is_classically_frustrated();
    if (is_signed)
      std::cerr << "WARNING: model has negative signs\n";
    if (is_classically_frustrated)
      std::cerr << "WARNING: model is classically frustrated\n";

    //
    // setup virtual lattice
    //

    bool is_bipartite = alps::set_parity(super_type::graph());
    vlat_.generate(rlat(), mp, mp.has_d_term());

    //
    // setup graph table
    //

    gtab_.clear();
    std::vector<double> weight;
    double rho = 0;
    site_iterator si, si_end;
    for (boost::tie(si, si_end) = alps::sites(rlat()); si != si_end; ++si) {
      looper::site_weight sw(mp.site(*si, rlat()));
      site_iterator vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) = virtual_sites(vlat(), rlat(), *si);
           vsi != vsi_end; ++vsi)
        for (int g = 1; g <= 3; ++g)
          if (alps::is_nonzero<1>(sw.v[g])) {
            gtab_.push_back(looper::site_graph(boost::get(
              looper::site_index_t(), vlat().graph(), *vsi), g));
            weight.push_back(sw.v[g]);
            rho += sw.v[g];
          }
    }
    bond_iterator bi, bi_end;
    for (boost::tie(bi, bi_end) = alps::bonds(rlat()); bi != bi_end; ++bi) {
      looper::bond_weight bw(mp.bond(*bi, rlat()));
      bond_iterator vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat(), rlat(), *bi);
           vbi != vbi_end; ++vbi)
        for (int g = 1; g <= 4; ++g)
          if (alps::is_nonzero<1>(bw.v[g])) {
            gtab_.push_back(looper::bond_graph(boost::get(
              looper::bond_index_t(), vlat().graph(), *vbi), g));
            weight.push_back(bw.v[g]);
            rho += bw.v[g];
          }
    }
    if (mp.has_d_term())
      for (boost::tie(si, si_end) = alps::sites(rlat()); si != si_end; ++si) {
        looper::bond_weight bw(mp.site(*si, rlat()));
        bond_iterator vbi, vbi_end;
        for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat(), rlat(), *si);
             vbi != vbi_end; ++vbi)
          for (int g = 1; g <= 4; ++g)
            if (alps::is_nonzero<1>(bw.v[g])) {
              gtab_.push_back(looper::bond_graph(boost::get(
                looper::bond_index_t(), vlat().graph(), *vbi), g));
              weight.push_back(bw.v[g]);
              rho += bw.v[g];
            }
      }
    r_graph.distribution().init(weight);

    //
    // init measurements
    //

    using alps::RealObservable;
    using alps::make_observable;

    if (is_signed) {
      super_type::measurements
        << RealObservable("Sign");
    }

    // unimproved super_type::measurements
    super_type::measurements
      << make_observable(
           RealObservable("Energy"), is_signed)
      << make_observable(
           RealObservable("Energy Density"), is_signed)
      << make_observable(
           RealObservable("Diagonal Energy Density"), is_signed)
      << make_observable(
           RealObservable("Energy Density^2"), is_signed)
      << make_observable(
           RealObservable("beta * Energy / sqrt(N)"), is_signed)
      << make_observable(
           RealObservable("beta * Energy^2"), is_signed)
      << make_observable(
           RealObservable("Susceptibility"), is_signed);
    if (is_bipartite)
      super_type::measurements
        << make_observable(
             RealObservable("Staggered Susceptibility"), is_signed);

    // improved measurements
    super_type::measurements
      << make_observable(
           RealObservable("Magnetization^2"), is_signed)
      << make_observable(
           RealObservable("Diagonal Energy Density (improved)"),
           is_signed);
    if (is_bipartite)
      super_type::measurements
        << make_observable(
             RealObservable("Staggered Magnetization^2"),
             is_signed);
    if (!is_classically_frustrated) {
      super_type::measurements
        << RealObservable("Uniform Generalized Magnetization^2")
        << RealObservable("Uniform Generalized Susceptibility");
      if (is_bipartite)
        super_type::measurements
          << RealObservable("Staggered Generalized Magnetization^2")
          << RealObservable("Staggered Generalized Susceptibility");
    }

  }

  virtual ~qmc_worker_base() {}

  virtual void dostep() { ++mcs_; }

  bool is_thermalized() const { return mcs_ >= mcs_therm_; }
  double work_done() const
  {
    return is_thermalized() ? (double(mcs_) / mcs_sweep_.min()) : 0.;
  }
  unsigned int mcs() const { return mcs_; }

  virtual void save(alps::ODump& od) const {
    super_type::save(od);
    od << mcs_;
  }
  virtual void load(alps::IDump& id) {
    super_type::load(id);
    id >> mcs_;
    if (super_type::where.empty()) super_type::measurements.compact();
  }

  const graph_type& rlat() const { return super_type::graph(); }
  const looper::virtual_lattice<graph_type>& vlat() const { return vlat_; }

private:
  unsigned int mcs_therm_;
  looper::integer_range<unsigned int> mcs_sweep_;
  looper::virtual_lattice<graph_type> vlat_;
  std::vector<looper::local_graph> gtab_;
  boost::variate_generator<alps::buffered_rng_base&, looper::random_choice<> >
    r_graph;

  // to be dumped/restored
  unsigned int mcs_;
};

struct path_integral {};
struct sse {};

template<class QMC,
  class MCRUN = alps::scheduler::LatticeModelMCRun<looper::graph_type> >
class qmc_worker;

template<class MCRUN>
class qmc_worker<looper::path_integral, MCRUN> : public qmc_worker_base<MCRUN>
{
public:
  typedef looper::path_integral qmc_type;
  typedef qmc_worker_base<MCRUN> super_type;
  typedef looper::local_operator<qmc_type> 
  
  class local_operator : public looper::local_operator_base
  {
  public:
    typedef looper::local_operator_base super_type;
    local_operator(bool is_bond, looper::operator_type type, unsigned int loc,
		   double time)
      : super_type(is_bond, type, loc), time_(time) {}
    double time() const { return time_; }
    void save(alps::ODump& dp) const { super_type::save(dp); dp << time_; }
    void load(alps::IDump& dp) { super_type::load(dp); dp >> time_; }
  private:
    double time_;
  };
  alps::ODump& operator<<(alps::ODump& dp, const local_operator& op)
  { op.save(dp); return dp; }
  alps::IDump& operator>>(alps::IDump& dp, local_operator& op)
  { op.load(dp); return dp; }

  qmc_worker(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n)
  {}

  void save(alps::ODump& od) const {
    super_type::save(od);
    od << operators << spins;
  }
  void load(alps::IDump& id) {
    super_type::load(id);
    id >> operators >> spins;
  }

private:
  std::vector<local_operator> operators;
  std::vector<int> spins;
  std::vector<local_operator> operators_prev; // no checkpointing
};

template<class MCRUN>
class qmc_worker<sse, MCRUN> : public qmc_worker_base<MCRUN>
{
public:
  typedef qmc_worker_base<MCRUN> super_type;

  qmc_worker(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n)
  {}

  void save(alps::ODump& od) const { super_type::save(od); }
  void load(alps::IDump& id) { super_type::load(id); }
};


template<class T>
inline void accumulate(const alps::ObservableSet& m_in, T& m_out)
{
  if (m_in.has("beta * Energy / sqrt(N)") && m_in.has("beta * Energy^2")) {
    alps::RealObsevaluator obse_e = m_in["beta * Energy / sqrt(N)"];
    alps::RealObsevaluator obse_e2 = m_in["beta * Energy^2"];
    alps::RealObsevaluator eval("Specific Heat");
    eval = (obse_e2 - obse_e * obse_e);
    m_out << eval;
  }
}

inline void accumulate(alps::scheduler::MCSimulation& sim)
{ accumulate(sim.get_measurements(), sim); }

#endif // LOOP_WORKER_H
