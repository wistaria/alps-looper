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

#include <looper/util.h>
#include <alps/alea.h>
#include <alps/scheduler.h>
#include <boost/regex.hpp>

namespace looper {

template<typename ESTIMATOR>
class evaluator : public alps::scheduler::MCRun
{
public:
  typedef ESTIMATOR              estimator_t;
  typedef alps::scheduler::MCRun super_type;

  evaluator(alps::ProcessList const& w, alps::Parameters const& p, int n)
    : super_type(w, p, n) {}

  void evaluate(alps::scheduler::MCSimulation& sim,
                alps::Parameters const& np) const
  {
    if (regex_match(sim.get_parameters().value_or_default("REPRESENTATION", ""),
                    boost::regex("QWL$"), boost::match_default)) {
      qwl_evaluate(sim.get_measurements(), sim.get_parameters(), np);
    } else {
      alps::ObservableSet m;
      evaluate(m, sim.get_measurements());
      for (alps::ObservableSet::const_iterator itr = m.begin(); itr != m.end();
           ++itr) sim.addObservable(*(itr->second));
    }
  }

  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in)
  { estimator_t::evaluate(m, m_in); }

  static void qwl_evaluate(alps::ObservableSet const& m,
                           alps::Parameters const& p,
                           alps::Parameters const& np)
  {
    if (!np.defined("T_MIN") && !p.defined("T_MIN"))
      boost::throw_exception(std::invalid_argument(
        "parameter T_MIN is not defined"));
    if (!np.defined("T_MAX") && !p.defined("T_MAX"))
      boost::throw_exception(std::invalid_argument(
        "parameter T_MAX is not defined"));
    if (!np.defined("DELTA_T") && !p.defined("DELTA_T"))
      boost::throw_exception(std::invalid_argument(
        "parameter DELTA_T is not defined"));
    double t_min = np.value_or_default("T_MIN", p.value_or_default("T_MIN", 0));
    double t_max = np.value_or_default("T_MAX", p.value_or_default("T_MAX", 0));
    double delta_t = np.value_or_default("DELTA_T",
                                         p.value_or_default("DELTA_T", 0));
  }
};

} // end namespace looper

#endif // LOOPER_EVALUATE_H
