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

#include "loop_worker.h"
#include <looper/weight.h>

qmc_worker_base::qmc_worker_base(alps::ProcessList const& w,
                                 alps::Parameters const& p, int n,
                                 bool is_path_integral)
  : super_type(w, p, n),
    mcs_sweep_(p.value_or_default("SWEEPS", "[65536:]")),
    mcs_therm_(p.value_or_default("THERMALIZATION", mcs_sweep_.min() >> 3)),
    mcs_(0),
    beta_(1.0 / static_cast<double>(p["T"])),
    chooser_(*engine_ptr)
{
  if (w == alps::ProcessList()) return;

  if (mcs_sweep_.min() < 1 || mcs_sweep_.min() > mcs_sweep_.max())
    boost::throw_exception(std::invalid_argument(
      "qmc_worker_base::qmc_worker_base() inconsistent MC steps"));

  if (beta_ < 0)
    boost::throw_exception(std::invalid_argument(
      "qmc_worker_base::qmc_worker_base() negative beta"));

  looper::model_parameter mp(p, *this);
  if (mp.is_signed()) std::cerr << "WARNING: model has negative signs\n";

  energy_offset_ = mp.energy_offset();

  vlat_.generate(graph(), mp, mp.has_d_term());
  if (mp.is_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";
  is_frustrated_ = mp.is_frustrated();
  double fs = p.value_or_default("FORCE_SCATTER", is_frustrated() ? 0.1 : 0);

  looper::weight_table wt(mp, rgraph(), vlattice(), fs);
  energy_offset_ += wt.energy_offset();
  chooser_.init(wt, is_path_integral);

  is_signed_ = mp.is_signed();
  if (mp.is_signed()) {
    bond_sign_.resize(num_bonds(vlat_));
    std::fill(bond_sign_.begin(), bond_sign_.end(), 0);
    looper::weight_table::bond_weight_iterator bi, bi_end;
    for (boost::tie(bi, bi_end) = wt.bond_weights(); bi != bi_end; ++bi)
      if (bi->second.sign < 0) bond_sign_[bi->first] = 1;
    site_sign_.resize(num_sites(vlat_));
    std::fill(site_sign_.begin(), site_sign_.end(), 0);
    looper::weight_table::site_weight_iterator si, si_end;
    for (boost::tie(si, si_end) = wt.site_weights(); si != si_end; ++si)
      if (si->second.sign < 0) site_sign_[si->first] = 1;
  }

  if (mp.has_field()) {
    field_.resize(0);
    site_iterator rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(rgraph()); rsi != rsi_end; ++rsi) {
      site_iterator vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) =
             virtual_sites(vlattice(), rgraph(), *rsi);
           vsi != vsi_end; ++vsi)
        field_.push_back(mp.site(*rsi, rgraph()).hz);
    }
  }

  use_improved_estimator_ =
    !(mp.has_field() || p.defined("DISABLE_IMPROVED_ESTIMATOR"));
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
