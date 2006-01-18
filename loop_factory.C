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
qmc_factory::make_task(const alps::ProcessList& w,
  const boost::filesystem::path& fn) const
{ return new alps::scheduler::MCSimulation(w, fn); }

alps::scheduler::MCSimulation*
qmc_factory::make_task(const alps::ProcessList& w,
  const boost::filesystem::path& fn, const alps::Parameters&) const
{ return new alps::scheduler::MCSimulation(w, fn); }

alps::scheduler::MCRun*
qmc_factory::make_worker(const alps::ProcessList& w,
  const alps::Parameters& p, int n) const
{
  if (p.defined("REPRESENTATION")) {
    map_type::const_iterator itr = creators_.find(p["REPRESENTATION"]);
    if (itr == creators_.end() || itr->second == 0)
      boost::throw_exception(std::runtime_error("unknown representation"));
    return itr->second->create(w, p, n);
  } else {
    if (creators_.size() == 1 && creators_.begin()->second)
      return creators_.begin()->second->create(w, p, n);
    else
      boost::throw_exception(std::runtime_error("no worker is registered"));
  }
  return 0;
}

void qmc_factory::print_copyright(std::ostream& os) const
{ looper::print_copyright(os); }

qmc_factory* qmc_factory::ptr_ = 0;
