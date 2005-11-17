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

std::ostream& operator<<(std::ostream& os, const looper::site_parameter& p)
{
  os << "C = " << p.c << ", Hx = " << p.hx << ", Hz = " << p.hz;
  return os;
}

std::ostream& operator<<(std::ostream& os, const looper::bond_parameter& p)
{
  os << "C = " << p.c << ", Jxy = " << p.jxy << ", Jz = " << p.jz;
  return os;
}

template<typename G>
void output(const alps::graph_helper<G>& gh,
            const looper::model_parameter& m)
{
  typedef G graph_type;
  typedef typename alps::graph_helper<G>::site_iterator site_iterator;
  typedef typename alps::graph_helper<G>::bond_iterator bond_iterator;

  // site parameters
  site_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = gh.sites(); vi != vi_end; ++vi) {
    std::cout << "site " << *vi << ": type = " << gh.site_type(*vi)
              << ", S = " << m.site(*vi, gh.graph()).s << std::endl;
  }

  // bond parameters
  bond_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = gh.bonds(); ei != ei_end; ++ei) {
    std::cout << "bond " << *ei << ": type = " << gh.bond_type(*ei)
              << ", " << m.bond(*ei, gh.graph()) << std::endl;
  }
}

qmc_worker_base::qmc_worker_base(const alps::ProcessList& w,
                                 const alps::Parameters& p, int n)
  : super_type(w, p, n),
    mcs_sweep_(p["SWEEPS"]),
    mcs_therm_(p.value_or_default("THERMALIZATION", mcs_sweep_.min() >> 3)),
    has_hz_(false),
    vlat_(), diag_graphs_(), offdiag_weights_(),
    r_graph_(*engine_ptr, looper::random_choice<>()),
    r_time_(*engine_ptr, boost::exponential_distribution<>()),
    mcs_(0)
{
  //
  // check parameters
  //

  if (mcs_sweep_.min() < mcs_therm_)
    boost::throw_exception(std::invalid_argument("too small SWEEPS"));

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

  bool is_bipartite = alps::set_parity(graph());
  vlat_.generate(rlat(), mp, mp.has_d_term());

  //
  // setup graph table and random number generators
  //

  output(*this, mp);
  looper::weight_table wt(mp, rlat(), vlat());
  wt.setup_graph_chooser(diag_graphs_, r_graph_, offdiag_weights_);
  r_time_.distribution() = boost::exponential_distribution<>(wt.rho());

  //
  // init measurements
  //

  using alps::RealObservable;
  using alps::make_observable;

  if (is_signed) {
    measurements << RealObservable("Sign");
  }

  // unimproved super_type::measurements
  measurements
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
    measurements
      << make_observable(
           RealObservable("Staggered Susceptibility"), is_signed);

  // improved measurements
  measurements
    << make_observable(
         RealObservable("Magnetization^2"), is_signed)
    << make_observable(
         RealObservable("Diagonal Energy Density (improved)"),
         is_signed);
  if (is_bipartite)
    measurements
      << make_observable(
           RealObservable("Staggered Magnetization^2"),
           is_signed);
  if (!is_classically_frustrated) {
    measurements
      << RealObservable("Uniform Generalized Magnetization^2")
      << RealObservable("Uniform Generalized Susceptibility");
    if (is_bipartite)
      measurements
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
  if (where.empty()) measurements.compact();
}
