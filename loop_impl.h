/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

// loop_impl.h - implementation of worker for QMC simulation

#ifndef LOOP_IMPL_H
#define LOOP_IMPL_H

#include <looper/copyright.h>
#include <looper/model.h>
#include <looper/path_integral.h>
#include <looper/sse.h>
#include <looper/measurement.h>
#include <alps/scheduler.h>

template<class QMC>
class qmc_worker
{
public:
  typedef QMC                      qmc;
  typedef typename qmc::graph_type graph_type;
  typedef typename qmc::model_type model_type;

  template<class G, class MDL>
  qmc_worker(const G& rg, const MDL& model, double beta, double fs,
             alps::ObservableSet& m) :
    param_(rg, model, beta, fs), config_()
  {
    using alps::RealObservable;
    using alps::make_observable;

    qmc::initialize(config_, param_);

// #ifndef NDEBUG
    if (is_signed()) std::cerr << "WARNING: model has negative signs\n";
    if (is_classically_frustrated())
      std::cerr << "WARNING: model is classically frustrated\n";
// #endif

    //
    // measurements
    //

    if (is_signed()) {
      m << RealObservable("Sign")
        << RealObservable("Sign (improved)");
    }

    // unimproved measurements
    m << make_observable(
           RealObservable("Energy"), is_signed())
      << make_observable(
           RealObservable("Energy Density"), is_signed())
      << make_observable(
           RealObservable("Diagonal Energy Density"), is_signed())
      << make_observable(
           RealObservable("Energy Density^2"), is_signed())
      << make_observable(
           RealObservable("beta * Energy / sqrt(N)"), is_signed())
      << make_observable(
           RealObservable("beta * Energy^2"), is_signed())
      << make_observable(
           RealObservable("Susceptibility"), is_signed());
    if (is_bipartite()) {
      m << make_observable(
             RealObservable("Staggered Susceptibility"), is_signed());
    }

    // improved measurements
    m << make_observable(
           RealObservable("Magnetization^2"),
           "Sign (improved)", double(), is_signed())
      << make_observable(
           RealObservable("Diagonal Energy Density (improved)"),
           "Sign (improved)", double(), is_signed());
    if (is_bipartite())
      m << make_observable(
             RealObservable("Staggered Magnetization^2"),
             "Sign (improved)", double(), is_signed());
    if (!is_classically_frustrated()) {
      m << RealObservable("Uniform Generalized Magnetization^2")
        << RealObservable("Uniform Generalized Susceptibility");
      if (is_bipartite())
        m << RealObservable("Staggered Generalized Magnetization^2")
          << RealObservable("Staggered Generalized Susceptibility");
    }
  }

  template<class RNG>
  void step(RNG& rng, alps::ObservableSet& m)
  {
    //
    // generate clusters
    //

    qmc::generate_loops(config_, param_, rng);

    //
    // measure improved quantities
    //

    if (is_signed()) {
      m["Sign (improved)"] << looper::sign_imp(config_, param_);
    }
    m["Diagonal Energy Density (improved)"] <<
      looper::energy_z_imp(config_, param_);
    m["Magnetization^2"] << looper::uniform_sz2_imp(config_, param_);
    if (is_bipartite())
      m["Staggered Magnetization^2"] <<
        looper::staggered_sz2_imp(config_, param_);

    if (!is_classically_frustrated()) {
      double gm2, gs;
      boost::tie(gm2, gs) =
        looper::uniform_generalized_susceptibility_imp(config_, param_);
      m["Uniform Generalized Magnetization^2"] << gm2;
      m["Uniform Generalized Susceptibility"] << gs;

      if (is_bipartite()) {
        double sgm2, sgs;
        boost::tie(sgm2, sgs) =
          looper::staggered_generalized_susceptibility_imp(config_, param_);
        m["Staggered Generalized Magnetization^2"] << sgm2;
        m["Staggered Generalized Susceptibility"] << sgs;
      }
    }

    //
    // flip clusters
    //

    qmc::flip_and_cleanup(config_, param_, rng);

    //
    // measure unimproved quantities
    //

    double sign = 1.0;
    if (is_signed()) {
      sign = looper::sign(config_, param_);
      m["Sign"] << sign;
    }
    
    double ez, exy, e2;
    boost::tie(ez, exy, e2) = looper::energy(config_, param_);
    m["Energy"] << param_.virtual_graph.num_real_vertices * (ez + exy);
    m["Energy Density"] << (ez + exy);
    m["Energy Density^2"] << e2;
    m["Diagonal Energy Density"] << ez;

    m["beta * Energy / sqrt(N)"] <<
      std::sqrt((double)param_.virtual_graph.num_real_vertices) *
      param_.beta * (ez + exy);
    m["beta * Energy^2"] <<
      param_.virtual_graph.num_real_vertices * looper::sqr(param_.beta) * e2;

    m["Susceptibility"] <<
      looper::uniform_susceptibility(config_, param_);
    if (is_bipartite()) {
      m["Staggered Susceptibility"] <<
        looper::staggered_susceptibility(config_, param_);
    }
  }

  void output_results(std::ostream& os, alps::ObservableSet& m)
  {
    using looper::mean; using looper::error;

    if (is_signed()) {
      os << mean(m["Sign"]) << ' ' << error(m["Sign"]) << ' '
         << mean(m["Sign (improved)"]) << ' '
         << error(m["Sign (improved)"]) << ' ';
    }
    os << mean(m["Energy"]) << ' ' << error(m["Energy"]) << ' '
       << mean(m["Energy Density"]) << ' ' << error(m["Energy Density"]) << ' '
       << mean(m["Specific Heat"]) << ' ' << error(m["Specific Heat"]) << ' '
       << mean(m["Magnetization^2"]) << ' '
       << error(m["Magnetization^2"]) << ' '
       << mean(m["Susceptibility"]) << ' '
       << error(m["Susceptibility"]) << ' ';
    if (is_bipartite()) {
      os << mean(m["Staggered Magnetization^2"]) << ' '
         << error(m["Staggered Magnetization^2"]) << ' '
         << mean(m["Staggered Susceptibility"]) << ' '
         << error(m["Staggered Susceptibility"]) << ' ';
    }
    if (!is_signed()) {
      os << mean(m["Uniform Generalized Magnetization^2"]) << ' '
         << error(m["Uniform Generalized Magnetization^2"]) << ' '
         << mean(m["Uniform Generalized Susceptibility"]) << ' '
         << error(m["Uniform Generalized Susceptibility"]) << ' ';
      if (is_bipartite()) {
        os << mean(m["Staggered Generalized Magnetization^2"]) << ' '
           << error(m["Staggered Generalized Magnetization^2"]) << ' '
           << mean(m["Staggered Generalized Susceptibility"]) << ' '
           << error(m["Staggered Generalized Susceptibility"]) << ' ';
      }
    }
  }

  void save(alps::ODump& od) const { config_.save(od); }
  void load(alps::IDump& id) { config_.load(id); }

protected:
  bool is_signed() const { return param_.model.is_signed(); }
  bool is_classically_frustrated() const
  { return param_.model.is_classically_frustrated(); }
  bool is_bipartite() const { return param_.is_bipartite; }

private:
  typename qmc::parameter_type param_;
  typename qmc::config_type    config_;
};

template<class T>
inline void accumulate(const alps::ObservableSet& m_in, T& m_out)
{
  alps::RealObsevaluator obse_e = m_in["beta * Energy / sqrt(N)"];
  alps::RealObsevaluator obse_e2 = m_in["beta * Energy^2"];
  alps::RealObsevaluator eval("Specific Heat");
  eval = (obse_e2 - obse_e * obse_e);
  m_out << eval;
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
    mdl_(p, graph(), operators(), model(), has_sign_problem()), mcs_(0),
    therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
    total_(therm_ + static_cast<unsigned int>(p["SWEEPS"])),
    strict_mcs_(p.defined("STRICT_MCS")),
    qmc_worker_(graph(), mdl_, 1.0 / static_cast<double>(p["T"]),
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

class factory : public alps::scheduler::Factory
{
  alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
    const boost::filesystem::path& fn) const
  {
    return new alps::scheduler::MCSimulation(w, fn);
  }
  alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
    const boost::filesystem::path& fn, const alps::Parameters&) const
  {
    return new alps::scheduler::MCSimulation(w, fn);
  }

  alps::scheduler::MCRun* make_worker(const alps::ProcessList& w,
                                      const alps::Parameters& p, int n) const
  {
    if (!p.defined("REPRESENTATION") ||
        p["REPRESENTATION"] == "path integral") {
      return new worker<qmc_worker<looper::path_integral<
        looper::virtual_graph<looper::parity_graph_type>,
        looper::model_parameter<> > > >(w, p, n);
    } else if (p["REPRESENTATION"] == "SSE") {
      return new worker<qmc_worker<looper::sse<
        looper::virtual_graph<looper::parity_graph_type>,
        looper::model_parameter<> > > >(w, p, n);
    } else {
      boost::throw_exception(std::invalid_argument("unknwon representation"));
    }
    return 0;
  }

  void print_copyright(std::ostream& os) const
  { looper::print_copyright(os); }
};

#endif // LOOP_IMPL_H
