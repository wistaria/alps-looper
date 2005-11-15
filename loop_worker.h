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
  typedef loop_config::graph_type                         graph_type;
  typedef alps::scheduler::LatticeModelMCRun<graph_type>  super_type;
  typedef looper::graph_traits<graph_type>::site_iterator site_iterator;
  typedef looper::graph_traits<graph_type>::bond_iterator bond_iterator;

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p, int n);
  virtual ~qmc_worker_base() {}

  virtual void dostep() { ++mcs_; }
  bool is_thermalized() const { return mcs_ >= mcs_therm_; }
  double work_done() const
  { return is_thermalized() ? (double(mcs_) / mcs_sweep_.min()) : 0.; }
  unsigned int mcs() const { return mcs_; }

  bool has_longitudinal_field() const { return has_hz_; }

  const graph_type& rlat() const { return super_type::graph(); }
  const looper::virtual_lattice<graph_type>& vlat() const { return vlat_; }
  unsigned int vsource(unsigned int b) const { return source(bond(b, vlat_)); }
  unsigned int vtarget(unsigned int b) const { return target(bond(b, vlat_)); }

  double advance() const { return r_time_(); }
  const looper::local_graph& choose_graph() const { return gtab_[r_graph_()]; }
  looper::local_graph choose_graph(const looper::location& loc) const
  {
    int g = (is_site(loc) || random() < otab_[pos(loc)]) ? 0 : 2;
    return looper::local_graph(g, loc);
  }

  virtual void save(alps::ODump& od) const;
  virtual void load(alps::IDump& id);

protected:
  double initialize(const looper::model_parameter& mp);

private:
  unsigned int mcs_therm_;
  looper::integer_range<unsigned int> mcs_sweep_;
  bool has_hz_;
  looper::virtual_lattice<graph_type> vlat_;
  std::vector<looper::local_graph> gtab_; // graph probability for diagonal conf
  std::vector<double> otab_; // graph probability for offdiagonal conf
  mutable boost::variate_generator<alps::buffered_rng_base&,
                                   looper::random_choice<> > r_graph_;
  mutable boost::variate_generator<alps::buffered_rng_base&,
                                   boost::exponential_distribution<> > r_time_;

  // to be dumped/restored
  unsigned int mcs_;
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
