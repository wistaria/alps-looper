/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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
    map_type::const_iterator itr = creators_.find(p["REPRESENTATION"]);
    if (itr == creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("unknown representation"));
    return itr->second->create(w, p, n);
  } else if (creators_.size() == 1 && creators_.begin()->second) {
    return creators_.begin()->second->create(w, p, n);
  } else {
    map_type::const_iterator itr = creators_.find("path integral");
    if (itr == creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("representation is not specified"));
    return itr->second->create(w, p, n);
  }
  return 0;
}

#ifdef HAVE_PARAPACK
alps::parapack::abstract_worker*
loop_factory::make_worker(alps::Parameters const& p, alps::ObservableSet& obs) const {
  if (p.defined("REPRESENTATION")) {
    map_type::const_iterator itr = creators_.find(p["REPRESENTATION"]);
    if (itr == creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("unknown representation"));
    return itr->second->create(p, obs);
  } else if (creators_.size() == 1 && creators_.begin()->second) {
    return creators_.begin()->second->create(p, obs);
  } else {
    map_type::const_iterator itr = creators_.find("path integral");
    if (itr == creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("representation is not specified"));
    return itr->second->create(p, obs);
  }
  return 0;
}
#endif // HAVE_PARAPACK

void loop_factory::print_copyright(std::ostream& os) const { looper::print_copyright(os); }

looper::abstract_evaluator* evaluator_factory::make_evaluator(const alps::Parameters& p) const {
  if (p.defined("REPRESENTATION")) {
    map_type::const_iterator itr = creators_.find(p["REPRESENTATION"]);
    if (itr == creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("unknown representation"));
    return itr->second->create();
  } else {
    if (creators_.size() == 1 && creators_.begin()->second)
      return creators_.begin()->second->create();
    else
      boost::throw_exception(std::runtime_error("representation is not specified"));
  }
  return 0;
}

void evaluator_factory::evaluate(alps::ObservableSet& obs, const alps::Parameters& p) const {
  looper::abstract_evaluator *eval = this->make_evaluator(p);
  eval->evaluate(obs, p, obs);
}

void evaluator_factory::print_copyright(std::ostream& os) const { looper::print_copyright(os); }

//
// initialization of static member pointer of factories
//

loop_factory* loop_factory::instance_ = 0;

evaluator_factory* evaluator_factory::instance_ = 0;
