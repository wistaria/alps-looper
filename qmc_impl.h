/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: qmc_impl.h 558 2003-11-12 14:02:13Z wistaria $
// qmc_impl.h - implementation of worker for QMC simulation

#ifndef QMC_IMPL_H
#define QMC_IMPL_H

#include <looper/copyright.h>
#include <looper/model.h>
#include <looper/path_integral.h>
#include <looper/sse.h>
#include <looper/measurement.h>
#include <looper/exact_diagonalization.h>

#include <alps/alea.h>
#include <alps/scheduler.h>

template<class QMC>
class qmc_worker
{
public:
  typedef QMC                      qmc;
  typedef typename qmc::graph_type graph_type;
  typedef alps::BasicSimpleObservable<double, alps::SimpleBinning<double> >
    measurement_type;
  typedef alps::RealObsevaluator   evaluator_type;

  template<class G, class MDL>
  qmc_worker(const G& rg, const MDL& model, double beta,
	     alps::ObservableSet& m) :
    param_(rg, model, beta), config_()
  {
    e_offset_ = looper::energy_offset(param_);

    qmc::initialize(config_, param_);

    // measurements
    m << measurement_type("diagonal energy");
    m << measurement_type("diagonal energy (improved)");
    m << measurement_type("off-diagonal energy");
    m << measurement_type("energy");
    m << measurement_type("energy^2");
    m << measurement_type("uniform magnetization^2");
    m << measurement_type("uniform magnetization^2 (improved)");
    m << measurement_type("uniform susceptibility");
    m << measurement_type("uniform susceptibility (improved)");
    m << measurement_type("staggered magnetization^2");
    m << measurement_type("staggered magnetization^2 (improved)");
    m << measurement_type("staggered susceptibility");
    m << measurement_type("staggered susceptibility (improved)");
    m << measurement_type("generalized magnetization^2 (improved)");
    m << measurement_type("generalized susceptibility (improved)");

    m << evaluator_type("specific heat");
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

    m.template
      get<measurement_type>("diagonal energy (improved)") <<
      looper::energy_z_imp(config_, param_);

    double m2 = looper::uniform_sz2_imp(config_, param_);
    m.template
      get<measurement_type>("uniform magnetization^2 (improved)") << m2;
    m.template
      get<measurement_type>("uniform susceptibility (improved)") <<
      param_.beta * m2;

    m.template
      get<measurement_type>("staggered magnetization^2 (improved)") << 
      looper::staggered_sz2_imp(config_, param_);

    double gm2, gs;
    boost::tie(gm2, gs) =
      looper::generalized_susceptibility_imp(config_, param_);
    m.template
      get<measurement_type>("generalized magnetization^2 (improved)") << gm2;
    m.template
      get<measurement_type>("generalized susceptibility (improved)") << gs;

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
    m.template get<measurement_type>("diagonal energy") << ez;
    m.template get<measurement_type>("off-diagonal energy") << exy;
    m.template get<measurement_type>("energy") << ez + exy;
    m.template get<measurement_type>("energy^2") << e2;

    m2 = param_.virtual_graph.num_real_vertices *
      looper::sqr(looper::uniform_sz(config_, param_));
    m.template get<measurement_type>("uniform magnetization^2") << m2;
    m.template get<measurement_type>("uniform susceptibility") <<
      param_.beta * m2;

    m.template get<measurement_type>("staggered magnetization^2") <<
      param_.virtual_graph.num_real_vertices *
      looper::sqr(looper::staggered_sz(config_, param_));
    m.template get<measurement_type>("staggered susceptibility") << 
      looper::staggered_susceptibility(config_, param_);
  }

  void accumulate(alps::ObservableSet& m)
  {
    evaluator_type obse_e = m.template get<measurement_type>("energy");
    evaluator_type obse_e2 = m.template get<measurement_type>("energy^2");
    m.template get<evaluator_type>("specific heat") = 
      (double)param_.virtual_graph.num_real_vertices *
      looper::sqr(param_.beta) * (obse_e2 - obse_e * obse_e);
  }

  static void output_results(std::ostream& os, alps::ObservableSet& m)
  {
    os << m.template get<measurement_type>("energy").mean() << ' '
       << m.template get<measurement_type>("energy").error() << ' '
       << m.template get<evaluator_type>("specific heat").mean() << ' '
       << m.template get<evaluator_type>("specific heat").error() << ' '
       << m.template get<measurement_type>("uniform susceptibility").mean() << ' '
       << m.template get<measurement_type>("uniform susceptibility").error() << ' '
       << m.template get<measurement_type>("staggered magnetization^2").mean() << ' '
       << m.template get<measurement_type>("staggered magnetization^2").error() << ' '
       << m.template get<measurement_type>("staggered susceptibility").mean() << ' '
       << m.template get<measurement_type>("staggered susceptibility").error() << ' ';
  }

  void save(alps::ODump& od) const { config_.save(od); }
  void load(alps::IDump& id) { config_.load(id); }

private:
  typename qmc::parameter_type param_;
  typename qmc::config_type    config_;
  bool                         is_bipartite;
  double                       e_offset_;
};


template<class QMC>
class worker : public alps::scheduler::LatticeModelMCRun<>
{
public:
  typedef alps::scheduler::LatticeModelMCRun<>::graph_type graph_type;

  worker(const alps::ProcessList& w, const alps::Parameters& p, int n) :
    alps::scheduler::LatticeModelMCRun<>(w, p, n),
    mdl_(p, graph(), simple_operators(), model()),
    mcs_(0), therm_(static_cast<unsigned int>(p["thermalization"])), 
    total_(static_cast<unsigned int>(p["MCS"])),
    qmc_(graph(), mdl_, 1.0 / static_cast<double>(p["temperature"]),
	 measurements) {}
  virtual ~worker() {}
    
  virtual void dostep() {
    ++mcs_;
    qmc_.step(random_01, measurements);
  }
  bool is_thermalized() const { return mcs_ >= therm_; }
  virtual double work_done() const {
    if (is_thermalized()) {
      return double(mcs_ - therm_) / (total_ - therm_);
    } else {
      return 0.;
    }
  }
  
  virtual void save(alps::ODump& od) const {
    od << mcs_;
    qmc_.save(od);
  }
  virtual void load(alps::IDump& id) {
    id >> mcs_;
    qmc_.load(id);
    if (where.empty()) measurements.compact();
  }

private:
  typename QMC::model_type mdl_;
  unsigned int mcs_, therm_, total_;
  qmc_worker<QMC> qmc_;
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
    if (!p.defined("representation") ||
	p["representation"] != "SSE") {
      return new worker<looper::path_integral<
                          looper::virtual_graph<looper::parity_graph_type>,
                          looper::xxz_model> >(w, p, n);
    } else {
      return new worker<looper::sse<
                          looper::virtual_graph<looper::parity_graph_type>,
	                  looper::xxz_model> >(w, p, n);
    }
  }

  void print_copyright(std::ostream& os) const
  { looper::print_copyright(os); }
};

#endif // QMC_IMPL_H
