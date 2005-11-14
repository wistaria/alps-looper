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
  : super_type(w, p, n),
    mcs_therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
    mcs_sweep_(p["SWEEPS"]), has_hz_(false),
    vlat_(), gtab_(), otab_(),
    r_graph_(*engine_ptr, looper::random_choice<>()),
    r_time_(*engine_ptr, boost::exponential_distribution<>()),
    mcs_(0)
{
  //
  // setup model
  //
  
  looper::model_parameter mp(p, *this);
  has_hz_ = mp.has_longitudinal_field();
  bool is_signed = mp.is_signed();
  bool is_classically_frustrated = mp.is_classically_frustrated();
  if (is_signed)
    std::cerr << "WARNING: model has negative signs\n";
  if (is_classically_frustrated)
    std::cerr << "WARNING: model is classically frustrated\n";
  
  //
  // setup virtual lattice
  //
  
  bool is_bipartite = alps::set_parity(super_type::graph());
  vlat_.generate(rlat(), mp, mp.has_d_term());
  
  //
  // setup graph table and random number generators
  //
  
  double rho = initialize(mp);
  r_time_.distribution() = boost::exponential_distribution<>(rho);
  
  //
  // init measurements
  //
  
  using alps::RealObservable;
  using alps::make_observable;
  
  if (is_signed) {
    super_type::measurements
      << RealObservable("Sign");
  }
  
  // unimproved super_type::measurements
  super_type::measurements
    << make_observable(
         RealObservable("Energy"), is_signed)
    << make_observable(
         RealObservable("Energy Density"), is_signed)
    << make_observable(
         RealObservable("Diagonal Energy Density"), is_signed)
    << make_observable(
         RealObservable("Energy Density^2"), is_signed)
    << make_observable(
         RealObservable("beta * Energy / sqrt(N)"), is_signed)
    << make_observable(
         RealObservable("beta * Energy^2"), is_signed)
    << make_observable(
         RealObservable("Susceptibility"), is_signed);
  if (is_bipartite)
    super_type::measurements
      << make_observable(
           RealObservable("Staggered Susceptibility"), is_signed);

  // improved measurements
  super_type::measurements
    << make_observable(
         RealObservable("Magnetization^2"), is_signed)
    << make_observable(
         RealObservable("Diagonal Energy Density (improved)"),
         is_signed);
  if (is_bipartite)
    super_type::measurements
      << make_observable(
           RealObservable("Staggered Magnetization^2"),
           is_signed);
  if (!is_classically_frustrated) {
    super_type::measurements
      << RealObservable("Uniform Generalized Magnetization^2")
      << RealObservable("Uniform Generalized Susceptibility");
    if (is_bipartite)
      super_type::measurements
        << RealObservable("Staggered Generalized Magnetization^2")
        << RealObservable("Staggered Generalized Susceptibility");
  }

}

void qmc_worker_base::save(alps::ODump& od) const {
  super_type::save(od);
  od << mcs_;
}

void qmc_worker_base::load(alps::IDump& id) {
  super_type::load(id);
  id >> mcs_;
  if (super_type::where.empty()) super_type::measurements.compact();
}

double qmc_worker_base::initialize(const looper::model_parameter& mp)
{
  gtab_.clear();
  otab_.clear();
  std::vector<double> weight;
  double rho = 0;
  site_iterator si, si_end;
  for (boost::tie(si, si_end) = sites(rlat()); si != si_end; ++si) {
    looper::site_weight sw(mp.site(*si, rlat()));
    site_iterator vsi, vsi_end;
    for (boost::tie(vsi, vsi_end) = virtual_sites(vlat_, rlat(), *si);
         vsi != vsi_end; ++vsi)
      for (int g = 0; g <= 2; ++g)
        if (alps::is_nonzero<1>(sw.v[g])) {
          gtab_.push_back(looper::site_graph(g, boost::get(
            looper::site_index_t(), vlat_.graph(), *vsi)));
          weight.push_back(sw.v[g]);
          rho += sw.v[g];
        }
  }
  bond_iterator bi, bi_end;
  for (boost::tie(bi, bi_end) = bonds(rlat()); bi != bi_end; ++bi) {
    looper::bond_weight bw(mp.bond(*bi, rlat()));
    bond_iterator vbi, vbi_end;
    for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat_, rlat(), *bi);
         vbi != vbi_end; ++vbi) {
      for (int g = 0; g <= 3; ++g)
        if (alps::is_nonzero<1>(bw.v[g])) {
          gtab_.push_back(looper::bond_graph(g, boost::get(
            looper::bond_index_t(), vlat_.graph(), *vbi)));
          weight.push_back(bw.v[g]);
          rho += bw.v[g];
        }
      if (alps::is_nonzero<1>(bw.v[0] + bw.v[2]))
        otab_.push_back(bw.v[0] / (bw.v[0] + bw.v[2]));
      else
        otab_.push_back(1);
    }
  }
  if (mp.has_d_term())
    for (boost::tie(si, si_end) = sites(rlat()); si != si_end; ++si) {
      looper::bond_weight bw(mp.site(*si, rlat()));
      bond_iterator vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat_, rlat(), *si);
           vbi != vbi_end; ++vbi)
        for (int g = 0; g <= 3; ++g)
          if (alps::is_nonzero<1>(bw.v[g])) {
            gtab_.push_back(looper::bond_graph(g, boost::get(
              looper::bond_index_t(), vlat_.graph(), *vbi)));
            weight.push_back(bw.v[g]);
            rho += bw.v[g];
          }
    }
  r_graph_.distribution().init(weight);
  return rho;
}
