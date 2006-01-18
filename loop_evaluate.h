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

#ifndef LOOP_EVALUATE_H
#define LOOP_EVALUATE_H

#include <looper/util.h>
#include <alps/alea.h>
#include <alps/scheduler.h>

template<typename ESTIMATOR>
class qmc_evaluator : public alps::scheduler::MCRun
{
public:
  typedef ESTIMATOR              estimator_t;
  typedef alps::scheduler::MCRun super_type;

  qmc_evaluator(alps::ProcessList const& w, alps::Parameters const& p, int n)
    : super_type(w, p, n) {}
  void load(alps::IDump& dp) { std::cerr << "loading information\n"; super_type::load(dp); dp >> info_; }

  void evaluate(alps::scheduler::MCSimulation& sim) const
  {
    alps::ObservableSet m;
    evaluate(m, info_, sim.get_measurements());
    for (alps::ObservableSet::const_iterator itr = m.begin(); itr != m.end();
         ++itr) sim << *(itr->second);
  }
  static void evaluate(alps::ObservableSet& m,
                       alps::Parameters const& info,
                       alps::ObservableSet const& m_in)
  {
    using looper::power2;
    if (info.defined("Number of Real Sites") &&
        info.defined("Inverse Temperature")) {
      int nrs = static_cast<int>(info["Number of Real Sites"]);
      double beta = static_cast<double>(info["Inverse Temperature"]);
      if (m_in.has("Energy") && m_in.has("Energy^2")) {
        alps::RealObsevaluator obse_e = m_in["Energy"];
        alps::RealObsevaluator obse_e2 = m_in["Energy^2"];
        alps::RealObsevaluator eval("Specific Heat");
        eval = power2(beta) * (obse_e2 - power2(obse_e)) / nrs;
        m << eval;
      }
      estimator_t::evaluate(m, beta, nrs, m_in);
    }
  }

private:
  alps::Parameters info_;
};

#endif // LOOP_EVALUATE_H
