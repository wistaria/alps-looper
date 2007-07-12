/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_config.h"
#include "loop_factory.h"
#include <looper/histogram.h>
#include <looper/montecarlo.h>
#include <alps/plot.h>
#include <boost/regex.hpp>

namespace {

class evaluator : public looper::abstract_evaluator {
public:
  typedef looper::measurement<loop_config::measurement_set>::type measurement_t;
  typedef measurement_t::evaluator evaluator_t;
  void pre_evaluate(alps::ObservableSet&, alps::Parameters const&,
    alps::ObservableSet const&) const {}
  void evaluate(alps::scheduler::MCSimulation& sim, alps::Parameters const&,
    boost::filesystem::path const&) const;
  void evaluate(alps::ObservableSet&, alps::Parameters const&,
    alps::ObservableSet const&) const {}
};

void evaluator::evaluate(alps::scheduler::MCSimulation& sim,
                         alps::Parameters const& np,
                         boost::filesystem::path const& f) const
{
  const alps::Parameters& p = sim.get_parameters();
  const alps::ObservableSet& m = sim.get_measurements();
  double nrs = alps::RealObsevaluator(m["Number of Sites"]).mean();

  typedef std::map<std::string, alps::plot::Set<double> > plot_set_type;
  plot_set_type plot_set;

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
    double maxg = g[0];
    for (int i = 1; i < g.size(); ++i) maxg = std::max(g[i] + lb*i, maxg);
    double z = 0;
    double sn = 0;
    double sn2 = 0;
    for (int i = 0; i < g.size(); ++i) {
      double c = std::exp(g[i] + lb*i - maxg);
      z += c;
      sn += i*c;
      sn2 += i*i*c;
    }
    double energy        = offset - t*sn/z;
    double free_energy   = offset - t * (std::log(z) + maxg);
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


//
// dynamic registration to the loop_factory
//

const bool pi_registered =
  loop_factory::instance()->register_evaluator<evaluator>
  ("path integral QWL");
const bool sse_registered =
  loop_factory::instance()->register_evaluator<evaluator>
  ("SSE QWL");

} // end namespace
