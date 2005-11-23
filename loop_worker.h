/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOP_WORKER_H
#define LOOP_WORKER_H

#include "loop_config.h"
#include <alps/alea.h>
#include <alps/scheduler.h>

class qmc_worker_base
  : public alps::scheduler::LatticeModelMCRun<loop_config::graph_type>
{
public:
  typedef alps::scheduler::LatticeModelMCRun<loop_config::graph_type>
    super_type;

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n),
      mcs_sweep_(p["SWEEPS"]),
      mcs_therm_(p.value_or_default("THERMALIZATION", mcs_sweep_.min() >> 3)),
      mcs_(0)
  {
    if (mcs_sweep_.min() < mcs_therm_)
      boost::throw_exception(std::invalid_argument(
        "qmc_worker_base::qmc_worker_base() too small SWEEPS"));
  }
  virtual ~qmc_worker_base() {}

  virtual void dostep() { ++mcs_; }
  bool is_thermalized() const { return mcs_ >= mcs_therm_; }
  double work_done() const
  { return is_thermalized() ? (double(mcs_) / mcs_sweep_.min()) : 0.; }
  unsigned int mcs() const { return mcs_; }

  const loop_config::graph_type& real_graph() const { return this->graph(); }
  loop_config::graph_type& real_graph() { return this->graph(); }

  virtual void save(alps::ODump& dp) const
  {
    super_type::save(dp);
    dp << mcs_;
  }
  virtual void load(alps::IDump& dp)
  {
    super_type::load(dp);
    dp >> mcs_;
    if (where.empty()) measurements.compact();
  }

private:
  looper::integer_range<unsigned int> mcs_sweep_;
  unsigned int mcs_therm_;
  unsigned int mcs_; // to be dumped/restored
};

template<class QMC> class qmc_worker;

template<class T>
inline void accumulate(const alps::ObservableSet& m_in, T& m_out)
{
  if (m_in.has("beta * Energy / sqrt(N)") && m_in.has("beta * Energy^2")) {
    alps::RealObsevaluator obse_e = m_in["beta * Energy / sqrt(N)"];
    alps::RealObsevaluator obse_e2 = m_in["beta * Energy^2"];
    alps::RealObsevaluator eval("Specific Heat");
    eval = (obse_e2 - obse_e * obse_e);
    m_out << eval;
  }
}

inline void accumulate(alps::scheduler::MCSimulation& sim)
{ accumulate(sim.get_measurements(), sim); }

#endif // LOOP_WORKER_H
