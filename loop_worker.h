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
#include <looper/model.h>
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

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p,
                  int n, bool is_path_integral = true);
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

  double beta() const { return beta_; }

  const lattice_graph_t& rgraph() const { return super_type::graph(); }
  const lattice_graph_t& vgraph() const { return vlat_.graph(); }
  const virtual_lattice& vlattice() const { return vlat_; }

  double energy_offset() const { return energy_offset_; }
  bool is_signed() const { return mp_.is_signed(); }
  bool is_classically_frustrated() const
  { return mp_.is_classically_frustrated(); }
  bool has_longitudinal_field() const { return mp_.has_longitudinal_field(); }
  bool has_d_term() const { return mp_.has_d_term(); }
  const looper::bond_parameter&
  rbond_parameter(const bond_descriptor& bd) const
  { return mp_.bond(bd, rgraph()); }
  const looper::site_parameter&
  rsite_parameter(const site_descriptor& sd) const
  { return mp_.site(sd, rgraph()); }

  loop_graph_t choose_graph() const { return chooser_.graph(); }
  loop_graph_t choose_diagonal(const location_t& loc, int c) const
  { return chooser_.diagonal(loc, c); }
  loop_graph_t choose_diagonal(const location_t& loc, int c0, int c1) const
  { return chooser_.diagonal(loc, c0, c1); }
  loop_graph_t choose_offdiagonal(const location_t& loc) const
  { return chooser_.offdiagonal(loc); }
  double advance() const { return chooser_.advance(); }
  double total_graph_weight() const { return chooser_.weight(); }

  void accumulate();
  void accumulate(const alps::ObservableSet& m_in, alps::ObservableSet& m_out);

  virtual void save(alps::ODump& dp) const;
  virtual void load(alps::IDump& dp);

private:
  looper::integer_range<unsigned int> mcs_sweep_;
  unsigned int mcs_therm_;
  unsigned int mcs_; // to be dumped/restored

  double beta_;

  looper::model_parameter mp_;
  double energy_offset_;
  virtual_lattice vlat_;
  graph_chooser chooser_;
};

template<class QMC> class qmc_worker;

#endif // LOOP_WORKER_H
