/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2003-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_EVALUATE_H
#define LOOPER_EVALUATE_H

#include "measurement.h"
#include "util.h"

#include <alps/alea.h>
#include <alps/plot.h>
#include <alps/scheduler.h>
#include <boost/regex.hpp>
#include <map>

namespace looper {

class abstract_evaluator : public alps::scheduler::MCRun
{
public:
  abstract_evaluator(alps::ProcessList const& w, alps::Parameters const& p,
                     int n) : alps::scheduler::MCRun(w, p, n) {}
  virtual ~abstract_evaluator() {}
  virtual void evaluate(alps::scheduler::MCSimulation& sim,
                        alps::Parameters const& np,
                        boost::filesystem::path const& f) const = 0;
};

template<typename ESTIMATOR>
class evaluator : public abstract_evaluator
{
public:
  typedef ESTIMATOR          estimator_t;
  typedef abstract_evaluator super_type;

  evaluator(alps::ProcessList const& w, alps::Parameters const& p, int n)
    : super_type(w, p, n) {}

  void evaluate(alps::scheduler::MCSimulation& sim,
                alps::Parameters const& np,
                boost::filesystem::path const& f) const
  {
    if (regex_match(sim.get_parameters().value_or_default("REPRESENTATION", ""),
                    boost::regex("(.*)QWL$"))) {
      qwl_evaluate(sim.get_measurements(), sim.get_parameters(), np, f);
    } else {
      alps::ObservableSet m;
      evaluate(m, sim.get_measurements());
      for (alps::ObservableSet::const_iterator itr = m.begin(); itr != m.end();
           ++itr) sim.addObservable(*(itr->second));
    }
  }

  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in)
  {
    energy_estimator::evaluate(m, m_in);
    estimator_t::evaluate(m, m_in);
  }

  static void qwl_evaluate(alps::ObservableSet const& m,
                           alps::Parameters const& p,
                           alps::Parameters const& np,
                           boost::filesystem::path const& f)
  {
    typedef std::map<std::string, alps::plot::Set<double> > plot_set_type;
    plot_set_type plot_set;

    double nrs = alps::RealObsevaluator(m["Number of Sites"]).mean();

    double t_min;
    if (np.defined("T_MIN"))
      t_min = np["T_MIN"];
    else if (p.defined("T_MIN"))
      t_min = p["T_MIN"];
    else
      boost::throw_exception(std::invalid_argument("T_MIN not defined"));

    double t_max;
    if (np.defined("T_MAX"))
      t_max = np["T_MAX"];
    else if (p.defined("T_MAX"))
      t_max = p["T_MAX"];
    else
      boost::throw_exception(std::invalid_argument("T_MAX not defined"));

    double t_delta;
    if (np.defined("T_DELTA"))
      t_delta = np["T_DELTA"];
    else if (p.defined("T_DELTA"))
      t_delta = p["T_DELTA"];
    else
      boost::throw_exception(std::invalid_argument("T_DELTA not defined"));

    // density of state
    std::valarray<double> g =
      alps::RealVectorObsevaluator(m["Partition Function Coefficient"]).mean();

    // energy
    double offset = alps::RealObsevaluator(m["Energy Offset"]).mean();
    plot_set["Energy Density"] = alps::plot::Set<double>();
    plot_set["Free Energy Density"] = alps::plot::Set<double>();
    plot_set["Entropy Density"] = alps::plot::Set<double>();
    plot_set["Specific Heat"] = alps::plot::Set<double>();

    for (double t = t_min; t < t_max + t_delta/2; t += t_delta) {
      double lb = -std::log(t);
      double maxx = g[0];
      for (int i = 1; i < g.size(); ++i) maxx = std::max(g[i] + lb*i, maxx);
      double z = 0;
      double sn = 0;
      double sn2 = 0;
      for (int i = 0; i < g.size(); ++i) {
        double c = std::exp(g[i] + lb*i - maxx);
        z += c;
        sn += i*c;
        sn2 += i*i*c;
      }
      double energy        = offset - t*sn/z;
      double free_energy   = offset - t * (std::log(z) + maxx);
      plot_set["Energy Density"] << t << energy / nrs;
      plot_set["Free Energy Density"] << t << free_energy / nrs;
      plot_set["Entropy Density"] << t << (energy - free_energy) / (t * nrs);
      plot_set["Specific Heat"] << t << (sn2/z - (sn/z) * (sn/z) - sn/z) / nrs;
    }

    // ouptut
    std::string prefix =
      regex_replace(f.string(), boost::regex("\\.out\\.xml$"), "");

    for (plot_set_type::const_iterator itr = plot_set.begin();
         itr != plot_set.end(); ++itr) {
      std::string long_name = itr->first;
      std::string short_name =
        regex_replace(long_name,
                      boost::regex("(\\sDensity$)|([A-Z])|(\\s)"),
                      "(?2\\l$2)(?3_)",
                      boost::match_default | boost::format_all);
      alps::plot::Plot<double> pl;
      pl.set_name(long_name + " versus Temperature (" + prefix + ")");
      pl.set_labels("Temperature", long_name);
      pl << itr->second;
      alps::oxstream oxs(prefix + ".plot." + short_name + ".xml");
      oxs << pl;
    }
  }
};

} // end namespace looper

#endif // LOOPER_EVALUATE_H
