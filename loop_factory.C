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

#include "loop_factory.h"
#include <looper/version.h>

alps::scheduler::MCSimulation*
loop_factory::make_task(const alps::ProcessList& w, const boost::filesystem::path& fn) const {
  return new alps::scheduler::MCSimulation(w, fn);
}

alps::scheduler::MCSimulation*
loop_factory::make_task(const alps::ProcessList& w, const boost::filesystem::path& fn,
  const alps::Parameters&) const {
  return new alps::scheduler::MCSimulation(w, fn);
}

alps::scheduler::MCRun*
loop_factory::make_worker(const alps::ProcessList& w, const alps::Parameters& p, int n) const {
  if (p.defined("REPRESENTATION")) {
    worker_map_type::const_iterator itr = worker_creators_.find(p["REPRESENTATION"]);
    if (itr == worker_creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("Unknown representation: " + p["REPRESENTATION"]));
    return itr->second->create(w, p, n);
  } else if (worker_creators_.size() == 1 && worker_creators_.begin()->second) {
    return worker_creators_.begin()->second->create(w, p, n);
  } else {
    worker_map_type::const_iterator itr = worker_creators_.find("path integral");
    if (itr == worker_creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("representation is not specified"));
    return itr->second->create(w, p, n);
  }
  return 0;
}

#ifdef HAVE_PARAPACK
alps::parapack::abstract_worker*
loop_factory::make_worker(alps::Parameters const& p, std::vector<alps::ObservableSet>& obs) const {
  if (p.defined("REPRESENTATION")) {
    worker_map_type::const_iterator itr = worker_creators_.find(p["REPRESENTATION"]);
    if (itr == worker_creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("Unknown representation: " + p["REPRESENTATION"]));
    return itr->second->create(p, obs);
  } else if (worker_creators_.size() == 1 && worker_creators_.begin()->second) {
    return worker_creators_.begin()->second->create(p, obs);
  } else {
    worker_map_type::const_iterator itr = worker_creators_.find("path integral");
    if (itr == worker_creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("representation is not specified"));
    return itr->second->create(p, obs);
  }
  return 0;
}
#endif // HAVE_PARAPACK

void loop_factory::print_copyright(std::ostream& os) const { looper::print_copyright(os); }

looper::abstract_evaluator* loop_factory::make_evaluator(const alps::Parameters& p) const {
  if (p.defined("REPRESENTATION")) {
    evaluator_map_type::const_iterator itr = evaluator_creators_.find(p["REPRESENTATION"]);
    if (itr == evaluator_creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("Unknown representation: " + p["REPRESENTATION"]));
    return itr->second->create();
  } else {
    if (evaluator_creators_.size() == 1 && evaluator_creators_.begin()->second)
      return evaluator_creators_.begin()->second->create();
    else
      boost::throw_exception(std::runtime_error("representation is not specified"));
  }
  return 0;
}

void loop_factory::pre_evaluate(std::vector<alps::ObservableSet>& obs,
  alps::Parameters const& p) const {
  looper::abstract_evaluator *eval = this->make_evaluator(p);
  BOOST_FOREACH(alps::ObservableSet& m, obs) eval->pre_evaluate(m, p, m);
}

void loop_factory::evaluate(std::vector<alps::ObservableSet>& obs,
  alps::Parameters const& p) const {
  looper::abstract_evaluator *eval = this->make_evaluator(p);
  BOOST_FOREACH(alps::ObservableSet& m, obs) eval->evaluate(m, p, m);
}

//
// initialization of static member pointer of factories
//

loop_factory* loop_factory::instance_ = 0;
