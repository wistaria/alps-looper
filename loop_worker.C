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
                                 const alps::Parameters& p, int n,
                                 bool is_path_integral)
  : super_type(w, p, n),
    mcs_sweep_(p.value_or_default("SWEEPS", "[65536:]")),
    mcs_therm_(p.value_or_default("THERMALIZATION", mcs_sweep_.min() >> 3)),
    mcs_(0),
    beta_(1.0 / static_cast<double>(p["T"])),
    chooser_(*engine_ptr)
{
  if (w != alps::ProcessList()) {
    if (mcs_sweep_.min() < 1 || mcs_sweep_.min() > mcs_sweep_.max())
      boost::throw_exception(std::invalid_argument(
        "qmc_worker_base::qmc_worker_base() inconsistent MC steps"));

    if (beta_ < 0)
      boost::throw_exception(std::invalid_argument(
        "qmc_worker_base::qmc_worker_base() negative beta"));

    looper::model_parameter mp(p, *this);
    if (mp.is_signed())
      std::cerr << "WARNING: model has negative signs\n";
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
      sign_.resize(0);
      looper::weight_table::bond_weight_iterator itr, itr_end;
      for (boost::tie(itr, itr_end) = wt.bond_weights(); itr != itr_end;
           ++itr)
        sign_.push_back(itr->second.sign);
    }

    has_field_ = mp.has_field();
    if (has_field()) {
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
  }
}

qmc_worker_base::~qmc_worker_base() {}

void qmc_worker_base::accumulate() { accumulate(measurements, measurements); }

void qmc_worker_base::accumulate(const alps::ObservableSet& m_in,
                                 alps::ObservableSet& m_out)
{
  using looper::sqr;
  if (m_in.has("Energy") && m_in.has("Energy^2")) {
    alps::RealObsevaluator obse_e = m_in["Energy"];
    alps::RealObsevaluator obse_e2 = m_in["Energy^2"];
    alps::RealObsevaluator eval("Specific Heat");
    eval = sqr(beta()) * (obse_e2 - sqr(obse_e)) / num_sites(rgraph());
    m_out << eval;
  }
}

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
