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

#include "loop_config.h"
#include "loop_factory.h"
#include <looper/cluster.h>
#include <looper/evaluator_impl.h>
#include <looper/model.h>
#include <looper/montecarlo.h>
#include <looper/permutation.h>
#include <alps/fixed_capacity_vector.h>
#include <boost/numeric/ublas/matrix.hpp>

namespace {

struct dummy_operator
{
  static bool is_offdiagonal() { return false; }
  static int pos() { return 0; }
  static bool is_site() { return false; }
};

class loop_worker
  : public alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
{
public:
  typedef looper::classical                        qmc_type;
  typedef alps::scheduler::LatticeModelMCRun<loop_config::lattice_graph_t>
                                                   super_type;
  typedef loop_config::lattice_graph_t             lattice_graph_t;
  typedef looper::virtual_lattice<lattice_graph_t> virtual_lattice;

  typedef std::vector<dummy_operator>   operator_string_t;
  typedef looper::union_find::node      cluster_fragment_t;
  typedef looper::cluster_info          cluster_info_t;

  typedef loop_config::estimator_t      estimator_t;
  typedef looper::measurement::estimate<estimator_t>::type estimate_t;

  loop_worker(alps::ProcessList const& w, alps::Parameters const& p, int n);
  void dostep();

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double work_done() const { return mcs.progress(); }

  void save(alps::ODump& dp) const
  { super_type::save(dp); dp << mcs << spins; }
  void load(alps::IDump& dp)
  { super_type::load(dp); dp >> mcs >> spins; }

  lattice_graph_t const& rgraph() const { return super_type::graph(); }
  lattice_graph_t const& vgraph() const { return vlattice.graph(); }

protected:
  void build();
  template<typename BIPARTITE, typename FIELD, typename IMPROVE>
  void flip();
  template<typename BIPARTITE, typename IMPROVE>
  void measure();

private:
  looper::mc_steps mcs;
  std::vector<int> spins;

  // parameters (no checkpointing)
  virtual_lattice vlattice;
  double beta;
  double energy_offset;
  bool use_improved_estimator;
  std::vector<double> coupling;
  boost::numeric::ublas::matrix<double,
    boost::numeric::ublas::column_major> prob;

  // working vectors (no checkpointing)
  std::vector<cluster_fragment_t> fragments;
  std::vector<cluster_info_t> clusters;
  std::vector<estimate_t> estimates;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::ProcessList const& w,
                         alps::Parameters const& p, int n)
  : super_type(w, p, n),
    beta(1.0 / alps::evaluate("T", p)),
    mcs(p)
{
  if (beta < 0)
    boost::throw_exception(std::invalid_argument(
      "loop_worker::loop_worker() negative temperature"));

  looper::model_parameter mp(p, *this);
  // if (!mp.is_classical())
  //   boost::throw_exception(std::invalid_argument(
  //     "loop_worker::loop_worker() not a classical Ising model"));
  if (mp.has_field())
    boost::throw_exception(std::invalid_argument(
      "longitudinal field is not supported yet"));
  if (mp.is_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";
  energy_offset = mp.energy_offset();

  vlattice.generate(rgraph(), mp, mp.has_d_term());

  coupling.resize(num_bonds(vgraph()));
  prob.resize(num_bonds(vgraph()), 2);
  int b = 0;
  bond_iterator rbi, rbi_end;
  for (boost::tie(rbi, rbi_end) = bonds(rgraph()); rbi != rbi_end; ++rbi) {
    double j = -mp.bond(*rbi, rgraph()).jz; // positive for ferromagnetic
    bond_iterator vbi, vbi_end;
    for (boost::tie(vbi, vbi_end) = virtual_bonds(vlattice, rgraph(), *rbi);
         vbi != vbi_end; ++vbi, ++b) {
      coupling[b] = j;
      prob(b, 0) = 1-std::exp(-2*beta*j); // for parallel
      prob(b, 1) = 1-std::exp( 2*beta*j); // for antiparallel
    }
  }

  // initialize configuration
  int nvs = num_sites(vgraph());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);

  // init measurements
  use_improved_estimator =
    !(mp.has_field() || p.defined("DISABLE_IMPROVED_ESTIMATOR"));
  measurements <<
    make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  measurements <<
    make_observable(alps::SimpleRealObservable("Number of Sites"));
  measurements <<
    make_observable(alps::SimpleRealObservable("Number of Clusters"));
  looper::energy_estimator::initialize(measurements, false);
  estimator_t::initialize(measurements, is_bipartite(), false,
                          use_improved_estimator);
}

void loop_worker::dostep()
{
  namespace mpl = boost::mpl;

  if (!mcs.can_work()) return;
  ++mcs;
  super_type::dostep();

  build();

  //   BIPARTITE    FIELD        IMPROVE
  flip<mpl::true_,  mpl::false_, mpl::true_ >();
  flip<mpl::true_,  mpl::false_, mpl::false_>();
  flip<mpl::false_, mpl::false_, mpl::true_ >();
  flip<mpl::false_, mpl::false_, mpl::false_>();

  //      BIPARTITE    IMPROVE
  measure<mpl::true_,  mpl::true_ >();
  measure<mpl::true_,  mpl::false_>();
  measure<mpl::false_, mpl::true_ >();
  measure<mpl::false_, mpl::false_>();
}


//
// cluster construction
//

void loop_worker::build()
{
  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(vgraph());
  fragments.resize(0); fragments.resize(nvs);

  int nvb = num_bonds(vgraph());
  for (int b = 0; b < nvb; ++b) {
    int s0 = vsource(b, vlattice);
    int s1 = vtarget(b, vlattice);
    if (random() < prob(b, spins[s0] ^ spins[s1])) unify(fragments, s0, s1);
  }

  // symmetrize spins
  alps::fixed_capacity_vector<int, loop_config::max_2s> r;
  site_iterator rsi, rsi_end;
  for (boost::tie(rsi, rsi_end) = sites(rgraph()); rsi != rsi_end; ++rsi) {
    site_iterator vsi, vsi_end;
    boost::tie(vsi, vsi_end) = virtual_sites(vlattice, rgraph(), *rsi);
    int offset = *vsi;
    int s2 = *vsi_end - *vsi;
    if (s2 > 1) {
      r.resize(s2);
      for (int i = 0; i < s2; ++i) r[i] = i;
      looper::restricted_random_shuffle(r.begin(), r.end(),
        spins.begin() + offset, spins.begin() + offset, random);
      for (int i = 0; i < s2; ++i) unify(fragments, i + offset, r[i] + offset);
    }
  }
}


//
// cluster flip
//

template<typename BIPARTITE, typename FIELD, typename IMPROVE>
void loop_worker::flip()
{
  if (!(is_bipartite() == BIPARTITE() &&
        false == FIELD() &&
        use_improved_estimator == IMPROVE())) return;

  int nvs = num_sites(vgraph());

  // assign cluster id
  int nc = 0;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) if (fi->is_root()) fi->id = nc++;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) fi->id = cluster_id(fragments, *fi);
  clusters.resize(0); clusters.resize(nc);

  if (IMPROVE()) estimates.resize(0); estimates.resize(nc);
  typename looper::measurement::accumulator<estimator_t, lattice_graph_t,
    time_t, cluster_fragment_t, BIPARTITE, IMPROVE>::type
    accum(estimates, fragments, vgraph());
  for (unsigned int s = 0; s < nvs; ++s) {
    accum.start(s, time_t(0), s, spins[s]);
    accum.term(s, time_t(1), s, spins[s]);
    accum.at_zero(s, s, spins[s]);
  }

  // determine whether clusters are flipped or not
  for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
       ci != clusters.end(); ++ci) ci->to_flip = (random() < 0.5);

  // flip spins
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id].to_flip) spins[s] ^= 1;

  // improved measurement
  if (IMPROVE()) {
    typename looper::measurement::collector<estimator_t, qmc_type, BIPARTITE,
      IMPROVE>::type coll;
    coll = std::accumulate(estimates.begin(), estimates.end(), coll);
    coll.commit(measurements, beta, num_sites(rgraph()), 0, 0);
  }
  measurements["Number of Clusters"] << (double)clusters.size();
}


//
// measurements
//

template<typename BIPARTITE, typename IMPROVE>
void loop_worker::measure()
{
  if (!(is_bipartite() == BIPARTITE() &&
        use_improved_estimator == IMPROVE())) return;

  int nrs = num_sites(rgraph());

  // energy
  double ene = energy_offset;
  int nvb = num_bonds(vgraph());
  for (int b = 0; b < nvb; ++b) {
    int s0 = vsource(b, vlattice);
    int s1 = vtarget(b, vlattice);
    ene -= 0.25 * coupling[b] * (0.5 - (spins[s0] ^ spins[s1]));
  }
  looper::energy_estimator::measure(measurements, beta, nrs, 0, 1, ene);

  looper::measurement::normal_estimator<estimator_t, qmc_type, BIPARTITE,
    IMPROVE>::type::measure(measurements, vgraph(), beta, nrs, 0, 1,
                            spins, operator_string_t(), spins);
}


//
// dynamic registration to the factories
//

const bool loop_registered =
  loop_factory::instance()->register_worker<loop_worker>("Ising");
const bool evaluator_registered = evaluator_factory::instance()->
  register_evaluator<looper::evaluator<loop_config::estimator_t> >("Ising");

} // end namespace
