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

#include "loop_worker.h"
#include <looper/weight.h>

qmc_worker_base::qmc_worker_base(const alps::ProcessList& w,
                                 const alps::Parameters& p, int n)
  : super_type(w, p, n),
    mcs_sweep_(p["SWEEPS"]),
    mcs_therm_(p.value_or_default("THERMALIZATION", mcs_sweep_.min() >> 3)),
    mcs_(0),
    mp_(p, *this), vlat_(graph(), mp_, mp_.has_d_term()),
    chooser_(*engine_ptr)
{
  if (mcs_sweep_.min() < mcs_therm_)
    boost::throw_exception(std::invalid_argument(
      "qmc_worker_base::qmc_worker_base() too small SWEEPS"));

  if (mp_.is_signed())
    std::cerr << "WARNING: model has negative signs\n";
  if (mp_.is_classically_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";

  energy_offset_ = 0;

  // energy offset in Hamiltonian
  for (site_iterator si = sites(rgraph()).first; si != sites(rgraph()).second;
       ++si) energy_offset_ += mp_.site(*si, rgraph()).c;
  for (bond_iterator bi = bonds(rgraph()).first; bi != bonds(rgraph()).second;
       ++bi) energy_offset_ += mp_.bond(*bi, rgraph()).c;
  if (mp_.has_d_term())
    for (site_iterator si = sites(rgraph()).first; si != sites(rgraph()).second;
         ++si) energy_offset_ += 0.5 * mp_.site(*si, rgraph()).s.get_twice();

  double fs = p.value_or_default("FORCE_SCATTER",
                                 mp_.is_classically_frustrated() ? 0.1 : 0);
  looper::weight_table wt(mp_, graph(), vlat_, fs);

  // energy offset in weight equation
  for (looper::weight_table::site_iterator si = wt.site_begin();
       si != wt.site_end(); ++si) energy_offset_ += si->second.offset;
  for (looper::weight_table::bond_iterator bi = wt.bond_begin();
       bi != wt.bond_end(); ++bi) energy_offset_ += bi->second.offset;

  chooser_.init(wt);
}

qmc_worker_base::~qmc_worker_base() {}

void qmc_worker_base::save(alps::ODump& dp) const
{
  super_type::save(dp);
  dp << mcs_;
}

void qmc_worker_base::load(alps::IDump& dp)
{
  super_type::load(dp);
  dp >> mcs_;
  if (where.empty()) measurements.compact();
}
