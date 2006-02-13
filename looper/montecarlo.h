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
  mc_steps operator++(int)
  { mc_steps tmp = *this; this->operator++(); return tmp; }

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

class wl_steps
{
public:
  wl_steps() {}
  wl_steps(alps::Parameters const& p)
    : steps_(0), logf_steps_(0), block_steps_(0), done_(false),
      num_logf_steps_(p.value_or_default("WANG_LANDAU_STEPS", 16)),
      num_block_steps_(p.value_or_default("BLOCK_STEPS", 10000))
  {}

  wl_steps& operator++()
  { ++steps_; ++logf_steps_; ++block_steps_; return *this; }
  wl_steps operator++(int)
  { wl_steps tmp = *this; this->operator++(); return tmp; }
  unsigned int operator()() const { return steps_; }
  bool can_work() const { true; }
  bool is_thermalized() const { return steps_; }
  double work_done() const { return done_ ? 1 : 0; }

  void save(alps::ODump& dp) const
  { dp << steps_ << logf_steps_ << block_steps_; }
  void load(alps::IDump& dp)
  { dp >> steps_ >> logf_steps_ >> block_steps_; }

private:
  unsigned int steps_;
  unsigned int logf_steps_;
  unsigned int block_steps_;
  bool done_;
  unsigned int num_logf_steps_;
  unsigned int num_block_steps_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::wl_steps const& mcs)
{ mcs.save(dp); return dp; }

inline alps::IDump& operator>>(alps::IDump& dp, looper::wl_steps& mcs)
{ mcs.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_MONTECARLO_H
