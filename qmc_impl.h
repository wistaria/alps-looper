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

// qmc_impl.h - implementation of worker for QMC simulation

#ifndef QMC_IMPL_H
#define QMC_IMPL_H

#include <looper/config.h>
#include <looper/copyright.h>
#include <looper/model.h>
#include <looper/path_integral.h>
#include <looper/sse.h>
#include <looper/measurement.h>

#ifdef HAVE_LAPACK
# include <looper/exact_diag.h>
#endif

#include <alps/alea.h>
#include <alps/scheduler.h>

template<class QMC>
class qmc_worker
{
public:
  BOOST_STATIC_CONSTANT(bool, is_qmc = true);

  typedef QMC                      qmc;
  typedef typename qmc::graph_type graph_type;
  typedef typename qmc::model_type model_type;
  typedef alps::RealObservable     measurement_type;
  typedef alps::RealObsevaluator   evaluator_type;

  template<class G, class MDL>
  qmc_worker(const G& rg, const MDL& model, double beta,
             alps::ObservableSet& m) :
    param_(rg, model, beta), config_()
  {
    e_offset_ = looper::energy_offset(param_);

    qmc::initialize(config_, param_);

    // measurements
    m << measurement_type("Energy");
    m << measurement_type("Energy Density");
    m << measurement_type("Energy Density^2");
    m << measurement_type("Magnetization^2");
    m << measurement_type("Susceptibility");
    m << measurement_type("Staggered Magnetization^2");
    m << measurement_type("Staggered Susceptibility");
    m << measurement_type("Generalized Magnetization^2");
    m << measurement_type("Generalized Susceptibility");

    m << evaluator_type("Specific Heat");
  }

  template<class RNG>
  void step(RNG& rng, alps::ObservableSet& m)
  {
    //
    // generate clusters
    //

    qmc::generate_loops(config_, param_, rng);

    //
    // measure improved quantities below
    //

    m.template get<measurement_type>("Magnetization^2") <<
      looper::uniform_sz2_imp(config_, param_);

    m.template get<measurement_type>("Staggered Magnetization^2") <<
      looper::staggered_sz2_imp(config_, param_);

    double gm2, gs;
    boost::tie(gm2, gs) =
      looper::generalized_susceptibility_imp(config_, param_);
    m.template get<measurement_type>("Generalized Magnetization^2") << gm2;
    m.template get<measurement_type>("Generalized Susceptibility") << gs;

    //
    // flip clusters
    //

    qmc::flip_and_cleanup(config_, param_, rng);

    //
    // measure unimproved quantities below
    //

    double ez, exy, e2;
    boost::tie(ez, exy, e2) = looper::energy(config_, param_);
    ez += e_offset_;
    m.template get<measurement_type>("Energy") <<
      param_.virtual_graph.num_real_vertices * (ez + exy);
    m.template get<measurement_type>("Energy Density") << ez + exy;
    m.template get<measurement_type>("Energy Density^2") << e2;

    double m2 = param_.virtual_graph.num_real_vertices *
      looper::sqr(looper::uniform_sz(config_, param_));
    m.template get<measurement_type>("Susceptibility") << param_.beta * m2;

    m.template get<measurement_type>("Staggered Susceptibility") <<
      looper::staggered_susceptibility(config_, param_);
  }

  void accumulate(alps::ObservableSet& m)
  {
    evaluator_type obse_e = m.template get<measurement_type>("Energy Density");
    evaluator_type obse_e2 =
      m.template get<measurement_type>("Energy Density^2");
    m.template get<evaluator_type>("Specific Heat") =
      (double)param_.virtual_graph.num_real_vertices *
      looper::sqr(param_.beta) * (obse_e2 - obse_e * obse_e);
  }

  static void output_results(std::ostream& os, alps::ObservableSet& m)
  {
    os << m.template get<measurement_type>("Energy").mean() << ' '
       << m.template get<measurement_type>("Energy").error() << ' '
       << m.template get<measurement_type>("Energy Density").mean() << ' '
       << m.template get<measurement_type>("Energy Density").error() << ' '
       << m.template get<evaluator_type>("Specific Heat").mean() << ' '
       << m.template get<evaluator_type>("Specific Heat").error() << ' '
       << m.template get<measurement_type>("Magnetization^2").mean() << ' '
       << m.template get<measurement_type>("Magnetization^2").error() << ' '
       << m.template get<measurement_type>("Susceptibility").mean() << ' '
       << m.template get<measurement_type>("Susceptibility").error() << ' '
       << m.template get<measurement_type>("Staggered Magnetization^2").mean() << ' '
       << m.template get<measurement_type>("Staggered Magnetization^2").error() << ' '
       << m.template get<measurement_type>("Staggered Susceptibility").mean() << ' '
       << m.template get<measurement_type>("Staggered Susceptibility").error() << ' '
       << m.template get<measurement_type>("Generalized Magnetization^2").mean() << ' '
       << m.template get<measurement_type>("Generalized Magnetization^2").error() << ' '
       << m.template get<measurement_type>("Generalized Susceptibility").mean() << ' '
       << m.template get<measurement_type>("Generalized Susceptibility").error() << ' ';
  }

  void save(alps::ODump& od) const { config_.save(od); }
  void load(alps::IDump& id) { config_.load(id); }

private:
  typename qmc::parameter_type param_;
  typename qmc::config_type    config_;
  double                       e_offset_;
};

#ifdef HAVE_LAPACK

template<class ED>
class ed_worker
{
public:
  BOOST_STATIC_CONSTANT(bool, is_qmc = false);

  typedef ED                      ed;
  typedef typename ed::graph_type graph_type;
  typedef typename ed::model_type model_type;
  typedef alps::BasicSimpleObservable<double, alps::NoBinning<double> >
    measurement_type;
  typedef alps::RealObsevaluator   evaluator_type;

  template<class G, class MDL>
  ed_worker(const G& rg, const MDL& model, double beta,
            alps::ObservableSet& m) :
    param_(rg, model, beta), config_(), done_(false)
  {
    // measurements
    m << measurement_type("Energy");
    m << measurement_type("Energy Density");
    m << measurement_type("Energy Density^2");
    m << measurement_type("Magnetization^2");
    m << measurement_type("Susceptibility");
    m << measurement_type("Staggered Magnetization^2");
    m << measurement_type("Staggered Susceptibility");

    m << evaluator_type("Specific Heat");
  }

  template<class RNG>
  void step(RNG& /* rng */, alps::ObservableSet& m)
  {
    if (done_) return;

    // generate Hamiltonian matrix
    ed::generate_matrix(param_, config_);

    // diagonalize matrix
    ed::diagonalize(param_, config_);

    // measure quantities
    double e, e2, c;
    boost::tie(e, e2, c) = ed::energy(param_, config_);
    double umag2, usus, smag2, ssus;
    boost::tie(umag2, usus, smag2, ssus) = ed::magnetization(param_, config_);

    m.template get<measurement_type>("Energy") <<
      (double)boost::num_vertices(param_.graph) * e;
    m.template get<measurement_type>("Energy Density") << e;
    m.template get<measurement_type>("Energy Density^2") << e2;
    m.template get<measurement_type>("Magnetization^2") << umag2;
    m.template get<measurement_type>("Susceptibility") << usus;
    m.template get<measurement_type>("Staggered Magnetization^2") << smag2;
    m.template get<measurement_type>("Staggered Susceptibility") << ssus;

    done_ = true;
  }

  void accumulate(alps::ObservableSet& m)
  {
    evaluator_type obse_e = m.template get<measurement_type>("Energy Density");
    evaluator_type obse_e2 =
      m.template get<measurement_type>("Energy Density^2");
    m.template get<evaluator_type>("Specific Heat") =
      (double)boost::num_vertices(param_.graph) *
      looper::sqr(param_.beta) * (obse_e2 - obse_e * obse_e);
  }

  static void output_results(std::ostream& os, alps::ObservableSet& m)
  {
    os << m.template get<measurement_type>("Energy").mean() << ' '
       << m.template get<measurement_type>("Energy").error() << ' '
       << m.template get<measurement_type>("Energy Density").mean() << ' '
       << m.template get<measurement_type>("Energy Density").error() << ' '
       << m.template get<evaluator_type>("Specific Heat").mean() << ' '
       << m.template get<evaluator_type>("Specific Heat").error() << ' '
       << m.template get<measurement_type>("Magnetization^2").mean() << ' '
       << m.template get<measurement_type>("Magnetization^2").error() << ' '
       << m.template get<measurement_type>("Susceptibility").mean() << ' '
       << m.template get<measurement_type>("Susceptibility").error() << ' '
       << m.template get<measurement_type>("Staggered Magnetization^2").mean() << ' '
       << m.template get<measurement_type>("Staggered Magnetization^2").error() << ' '
       << m.template get<measurement_type>("Staggered Susceptibility").mean() << ' '
       << m.template get<measurement_type>("Staggered Susceptibility").error() << ' ';
  }

  void save(alps::ODump& /* od */) const {}
  void load(alps::IDump& /* id */) {}

private:
  typename ed::parameter_type param_;
  typename ed::config_type    config_;
  bool                        done_;
};

#endif // HAVE_LAPACK

template<class QMC_WORKER>
class worker : public alps::scheduler::LatticeModelMCRun<>
{
public:
  typedef alps::scheduler::LatticeModelMCRun<>::graph_type graph_type;

  worker(const alps::ProcessList& w, const alps::Parameters& p, int n) :
    alps::scheduler::LatticeModelMCRun<>(w, p, n),
    mdl_(p, graph(), simple_operators(), model()),
    mcs_(0), therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
    total_(therm_ + static_cast<unsigned int>(p["SWEEPS"])),
    qmc_worker_(graph(), mdl_, 1.0 / static_cast<double>(p["T"]),
                measurements) {}
  virtual ~worker() {}

  virtual void dostep() {
    if (QMC_WORKER::is_qmc || is_thermalized())
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
                          looper::xxz_model> > >(w, p, n);
    } else if (p["REPRESENTATION"] == "SSE") {
      return new worker<qmc_worker<looper::sse<
                          looper::virtual_graph<looper::parity_graph_type>,
                          looper::xxz_model> > >(w, p, n);
#ifdef HAVE_LAPACK
    } else if (p["REPRESENTATION"] == "exact diagonalization") {
      return new worker<ed_worker<looper::exact_diagonalization<
                          looper::parity_graph_type,
                          looper::xxz_model> > >(w, p, n);
#endif // HAVE_LAPACK
    } else {
      boost::throw_exception(std::invalid_argument("unknwon representation"));
    }
    return 0;
  }

  void print_copyright(std::ostream& os) const
  { looper::print_copyright(os); }
};

#endif // QMC_IMPL_H
