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

#ifndef LOOPER_EVALUATOR_IMPL_H
#define LOOPER_EVALUATOR_IMPL_H

#include "evaluator.h"
#include "measurement.h"

namespace looper {

template<typename ESTIMATOR>
class evaluator : public abstract_evaluator
{
public:
  typedef ESTIMATOR estimator_t;
  void evaluate(alps::scheduler::MCSimulation& sim, alps::Parameters const&,
                boost::filesystem::path const&) const
  {
    alps::ObservableSet m;
    evaluate(m, sim.get_measurements());
    for (alps::ObservableSet::const_iterator itr = m.begin(); itr != m.end();
         ++itr) sim.addObservable(*(itr->second));
  }

  void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in) const
  {
    looper::energy_estimator::evaluate(m, m_in);
    estimator_t::evaluate(m, m_in);
  }
};

} // end namespace looper

#endif // LOOPER_EVALUATOR_IMPL_H
