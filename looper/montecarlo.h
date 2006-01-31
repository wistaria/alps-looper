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

#ifndef LOOPER_MONTECARLO_H
#define LOOPER_MONTECARLO_H

#include "util.h"
#include <alps/osiris.h>
#include <alps/scheduler.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class mc_steps
{
public:
  mc_steps() : mcs_(0) {}
  mc_steps(alps::Parameters const& p)
    : mcs_(0),
      sweep_(p.value_or_default("SWEEPS", "[65536:]")),
      therm_(p.value_or_default("THERMALIZATION", sweep_.min() >> 3))
  {
    if (sweep_.min() < 1 || sweep_.min() > sweep_.max())
      boost::throw_exception(std::invalid_argument(
        "qmc_worker::qmc_worker() inconsistent MC steps"));
  }

  mc_steps& operator++() { ++mcs_; return *this; }
  mc_steps operator++(int) { mc_steps tmp = *this; ++mcs_; return tmp; }

  unsigned int operator()() const { return mcs_; }
  bool can_work() const
  { return mcs_ < therm_ || mcs_ - therm_ < sweep_.max(); }
  bool is_thermalized() const { return mcs_ >= therm_; }
  double work_done() const
  { return is_thermalized() ? (double(mcs_ - therm_) / sweep_.min()) : 0; }

  unsigned int num_thermalization() const { return therm_; }
  std::pair<unsigned int, unsigned int> num_sweeps() const
  { return std::make_pair(sweep_.min(), sweep_.max()); }

  void save(alps::ODump& dp) const { dp << mcs_ << sweep_ << therm_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> sweep_ >> therm_; }

private:
  unsigned int mcs_;
  looper::integer_range<unsigned int> sweep_;
  unsigned int therm_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::mc_steps const& mcs)
{ mcs.save(dp); return dp; }

inline alps::IDump& operator>>(alps::IDump& dp, looper::mc_steps& mcs)
{ mcs.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

template<typename GRAPH>
class mc_worker : public alps::scheduler::LatticeModelMCRun<GRAPH>
{
public:
  typedef GRAPH                                               lattice_graph_t;
  typedef alps::scheduler::LatticeModelMCRun<lattice_graph_t> super_type;

  mc_worker(alps::ProcessList const& w, alps::Parameters const& p, int n)
    : super_type(w, p, n), mcs_(p) {}
  virtual ~mc_worker() {}

  virtual void dostep() { ++mcs_; }
  virtual void evaluate() = 0;

  virtual void save(alps::ODump& dp) const { super_type::save(dp); dp << mcs_; }
  virtual void load(alps::IDump& dp) { super_type::load(dp); dp >> mcs_; }

  bool can_work() const { return mcs_.can_work(); }
  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double work_done() const { return mcs_.work_done(); }
  unsigned int mcs() const { return mcs_(); }

  unsigned int num_thermalization() const { return mcs_.num_thermalization(); }
  std::pair<unsigned int, unsigned int> num_sweeps() const
  { return mcs_.num_sweeps(); }

private:
  looper::mc_steps mcs_;
};

} // end namepspace looper

#endif // LOOPER_MONTECARLO_H
