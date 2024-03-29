/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2003-2010 by Synge Todo <wistaria@comp-phys.org>
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
#include <boost/foreach.hpp>

namespace looper {

template<typename MEASUREMENT_SET>
class evaluator : public abstract_evaluator {
public:
  typedef typename measurement<MEASUREMENT_SET>::type measurement_t;
  evaluator() {}
  evaluator(alps::Parameters const& p) : params(p) {}
  void pre_evaluate(alps::ObservableSet& m, alps::Parameters const& p,
    alps::ObservableSet const& m_in) const {
    measurement_t::pre_evaluator::pre_evaluate(m, p, m_in);
  }
  void evaluate(alps::scheduler::MCSimulation& sim, alps::Parameters const& p,
    boost::filesystem::path const&) const {
    alps::ObservableSet m;
    measurement_t::evaluator::evaluate(m, p, sim.get_measurements());
    BOOST_FOREACH(alps::ObservableSet::iterator::value_type const& v, m)
      sim.addObservable(*(v.second));
  }
  void evaluate(alps::ObservableSet& m, alps::Parameters const& p,
    alps::ObservableSet const& m_in) const {
    measurement_t::evaluator::evaluate(m, p, m_in);
  }
  void evaluate(alps::ObservableSet& m) const {
    measurement_t::evaluator::evaluate(m, params, m);
  }
private:
  alps::Parameters params;
};

} // end namespace looper

#endif // LOOPER_EVALUATOR_IMPL_H
