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

qmc_worker_base::qmc_worker_base(const alps::ProcessList& w,
                                 const alps::Parameters& p, int n)
  : lattice_model_mcrun(w, p, n),
    virtual_lattice_adaptor(lattice_model_mcrun::graph()),
    mcs_sweep_(p["SWEEPS"]),
    mcs_therm_(p.value_or_default("THERMALIZATION", mcs_sweep_.min() >> 3)),
    mcs_(0),
    mp(p, *this), chooser(*engine_ptr)
{
  if (mcs_sweep_.min() < mcs_therm_)
    boost::throw_exception(std::invalid_argument(
      "qmc_worker_base::qmc_worker_base() too small SWEEPS"));

  mp.set_parameters(p, *this);
  if (mp.is_signed())
    std::cerr << "WARNING: model has negative signs\n";
  if (mp.is_classically_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";

  energy_offset = 0;
  site_iterator si, si_end;
  for (boost::tie(si, si_end) = sites(rgraph()); si != si_end; ++si)
    energy_offset += mp.site(*si, rgraph()).c;
  bond_iterator bi, bi_end;
  for (boost::tie(bi, bi_end) = bonds(rgraph()); bi != bi_end; ++bi)
    energy_offset += mp.bond(*bi, rgraph()).c;
  if (mp.has_d_term())
    for (boost::tie(si, si_end) = sites(rgraph()); si != si_end; ++si)
      energy_offset += 0.5 * mp.site(*si, rgraph()).s.get_twice();

  virtual_lattice_adaptor::init(mp);

  chooser.init(looper::weight_table(mp, rgraph(), vlattice()));
}

qmc_worker_base::~qmc_worker_base() {}

void qmc_worker_base::save(alps::ODump& dp) const
{
  lattice_model_mcrun::save(dp);
  dp << mcs_;
}

void qmc_worker_base::load(alps::IDump& dp)
{
  lattice_model_mcrun::load(dp);
  dp >> mcs_;
  if (where.empty()) measurements.compact();
}
