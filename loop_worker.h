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
  : public alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
{
public:
  typedef alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
    super_type;
  typedef loop_config::lattice_graph_t             lattice_graph_t;
  typedef loop_config::loop_graph_t            loop_graph_t;
  typedef loop_graph_t::location_t             location_t;
  typedef looper::virtual_lattice<lattice_graph_t> virtual_lattice;
  typedef looper::graph_chooser<loop_graph_t, super_type::engine_type>
    graph_chooser;

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p, int n);
  virtual ~qmc_worker_base();

  virtual void dostep() { ++mcs_; }
  bool can_work() const { return mcs_ < mcs_therm_ + mcs_sweep_.max(); }
  bool is_thermalized() const { return mcs_ >= mcs_therm_; }
  double work_done() const
  {
    return is_thermalized() ?
      (double(mcs_ - mcs_therm_) / mcs_sweep_.min()) : 0.;
  }
  unsigned int mcs() const { return mcs_; }

  const lattice_graph_t& rgraph() const { return super_type::graph(); }
  const lattice_graph_t& vgraph() const { return vlat_.graph(); }
  const virtual_lattice& vlattice() const { return vlat_; }

  loop_graph_t choose_diagonal() const { return chooser.diagonal(); }
  loop_graph_t choose_offdiagonal(const location_t& loc) const
  { return chooser.offdiagonal(loc); }
  double advance() const { return chooser.advance(); }

  virtual void save(alps::ODump& dp) const;
  virtual void load(alps::IDump& dp);

private:
  looper::integer_range<unsigned int> mcs_sweep_;
  unsigned int mcs_therm_;
  unsigned int mcs_; // to be dumped/restored

protected:
  looper::model_parameter mp;
  double energy_offset;

private:
  virtual_lattice vlat_;

protected:
  graph_chooser chooser;
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
