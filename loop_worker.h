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

#include <alps/alea.h>
#include <alps/scheduler.h>
#include <looper.h>

template<class MCRUN = alps::scheduler::LatticeModelMCRun<looper::graph_type> >
class qmc_worker_base : public MCRUN
{
public:
  typedef MCRUN                           super_type;
  typedef typename super_type::graph_type graph_type;

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n),
      mcs_therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
      mcs_sweep_(p["SWEEPS"]),
      vlat_()
  {
    using alps::RealObservable;
    using alps::make_observable;

    // set virtual graph and weights
    looper::model_parameter mp(p, *this);

    // warning
    if (super_type::is_signed())
      std::cerr << "WARNING: model has negative signs\n";
    if (super_type::is_classically_frustrated())
      std::cerr << "WARNING: model is classically frustrated\n";

    //
    // init measurements
    //
    if (super_type::is_signed()) {
      super_type::measurements
	<< RealObservable("Sign");
    }

    // unimproved super_type::measurements
    super_type::measurements
      << make_observable(
           RealObservable("Energy"), super_type::is_signed())
      << make_observable(
           RealObservable("Energy Density"), super_type::is_signed())
      << make_observable(
           RealObservable("Diagonal Energy Density"), super_type::is_signed())
      << make_observable(
           RealObservable("Energy Density^2"), super_type::is_signed())
      << make_observable(
           RealObservable("beta * Energy / sqrt(N)"), super_type::is_signed())
      << make_observable(
           RealObservable("beta * Energy^2"), super_type::is_signed())
      << make_observable(
           RealObservable("Susceptibility"), super_type::is_signed());
    if (super_type::is_bipartite())
      super_type::measurements
	<< make_observable(
             RealObservable("Staggered Susceptibility"), super_type::is_signed());

    // improved measurements
    super_type::measurements
      << make_observable(
           RealObservable("Magnetization^2"), super_type::is_signed())
      << make_observable(
           RealObservable("Diagonal Energy Density (improved)"),
           super_type::is_signed());
    if (super_type::is_bipartite())
      super_type::measurements
	<< make_observable(
             RealObservable("Staggered Magnetization^2"),
             super_type::is_signed());
    if (!super_type::is_classically_frustrated()) {
      super_type::measurements
	<< RealObservable("Uniform Generalized Magnetization^2")
        << RealObservable("Uniform Generalized Susceptibility");
      if (super_type::is_bipartite())
        super_type::measurements
	  << RealObservable("Staggered Generalized Magnetization^2")
          << RealObservable("Staggered Generalized Susceptibility");
    }

  }

  virtual ~qmc_worker_base() {}

  void dostep() { ++mcs_; }

  bool is_thermalized() const { return mcs_ >= mcs_therm_; }
  double work_done() const
  {
    return is_thermalized() ? (double(mcs_) / mcs_sweep_.min()) : 0.;
  }

  void save(alps::ODump& od) const {
    od << mcs_;
    super_type::save(od);
  }
  virtual void load(alps::IDump& id) {
    id >> mcs_;
    super_type::load(id);
    if (super_type::where.empty()) super_type::measurements.compact();
  }

  const looper::virtual_lattice<graph_type>& vlattice() const { return vlat_; }

private:
  unsigned int mcs_therm_;
  looper::integer_range<unsigned int> mcs_sweep_;
  looper::virtual_lattice<graph_type> vlat_;

  // to be dumped/restored
  unsigned int mcs_;
};

struct path_integral {};
struct sse {};

template<class QMC,
  MCRUN = alps::scheduler::LatticeModelMCRun<looper::graph_type> >
class qmc_worker;

template<class MCRUN>
class qmc_worker<path_integral, MCRUN> : public qmc_worker_base<MCRUN>
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

template<class QMC_WORKER>
class worker : public alps::scheduler::LatticeModelMCRun<>
{
public:
  typedef alps::scheduler::LatticeModelMCRun<>::graph_type graph_type;

  worker(const alps::ProcessList& w, const alps::Parameters& p, int n) :
    alps::scheduler::LatticeModelMCRun<>(w, p, n),
    mdl_(p, this->graph(), this->inhomogeneous_sites(),
         this->inhomogeneous_bonds(), *this, has_sign_problem()), mcs_(0),
    therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
    total_(therm_ + static_cast<unsigned int>(p["SWEEPS"])),
    strict_mcs_(p.defined("STRICT_MCS")),
    qmc_worker_(this->graph(), mdl_, 1.0 / static_cast<double>(p["T"]),
                p.value_or_default("FORCE_SCATTER", 0.0),
                measurements)
  { if (p.defined("FIXED_SEED")) random.seed(parms["FIXED_SEED"]); }
  virtual ~worker() {}

  virtual void dostep() {
    if (!strict_mcs_ || mcs_ < total_)
      qmc_worker_.step(random_01, measurements);
    ++mcs_;
  }
  bool is_thermalized() const { return mcs_ >= therm_; }
  virtual double work_done() const {
    return is_thermalized() ?
      double(mcs_ - therm_) / (total_ - therm_) : 0.0;
  }

  virtual void save(alps::ODump& od) const {
    od << mcs_;
    qmc_worker_.save(od);
  }
  virtual void load(alps::IDump& id) {
    id >> mcs_;
    qmc_worker_.load(id);
    if (where.empty()) measurements.compact();
  }

private:
  typename QMC_WORKER::model_type mdl_;
  unsigned int mcs_, therm_, total_;
  bool strict_mcs_;
  QMC_WORKER qmc_worker_;
};
