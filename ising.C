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
#include <looper/temperature.h>
#include <boost/numeric/ublas/matrix.hpp>

namespace {

struct dummy_operator
{
  static bool is_offdiagonal() { return false; }
  static int pos() { return 0; }
  static bool is_site() { return false; }
  static bool is_bond() { return false; }
};

class loop_worker : public alps::scheduler::MCRun
{
public:
  typedef looper::classical             mc_type;
  typedef alps::scheduler::MCRun        super_type;

  typedef looper::lattice_helper<loop_config::lattice_graph_t>
                                        lattice_t;
  typedef std::vector<dummy_operator>   operator_string_t;

  typedef looper::union_find::node      cluster_fragment_t;
  typedef looper::cluster_info          cluster_info_t;

  typedef looper::estimator<loop_config::measurement_set, mc_type, lattice_t,
    time_t>::type                       estimator_t;

  loop_worker(alps::ProcessList const& w, alps::Parameters const& p, int n);
  void dostep();

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double work_done() const { return mcs.progress(); }

  void save(alps::ODump& dp) const
  { super_type::save(dp); dp << mcs << spins; }
  void load(alps::IDump& dp)
  { super_type::load(dp); dp >> mcs >> spins; }

protected:
  void build();
  template<typename FIELD, typename IMPROVE>
  void flip();
  void measure();

private:
  // helpers
  lattice_t lattice;
  looper::spinmodel_helper model;

  // parameters
  looper::temperature temperature;
  double beta;
  double energy_offset;
  bool use_improved_estimator;
  std::vector<double> field;
  std::vector<double> coupling;
  boost::numeric::ublas::matrix<double> prob;

  // configuration (checkpoint)
  looper::mc_steps mcs;
  std::vector<int> spins;

  // observables
  alps::ObservableSet& obs;
  estimator_t estimator;

  // working vectors
  std::vector<cluster_fragment_t> fragments;
  std::vector<cluster_info_t> clusters;
  std::vector<looper::estimate<estimator_t>::type> estimates;
  std::vector<int> perm;
};


//
// member functions of loop_worker
//

loop_worker::loop_worker(alps::ProcessList const& w,
                         alps::Parameters const& p, int n)
  : super_type(w, p, n), lattice(p), model(p, lattice),
    temperature(p), mcs(p), obs(measurements)
{
  if (temperature.annealing_steps() > mcs.thermalization())
    boost::throw_exception(std::invalid_argument(
      "annealing steps are longer than thermalization steps"));

  energy_offset = model.energy_offset();
  use_improved_estimator =
    !model.has_field() && !p.defined("DISABLE_IMPROVED_ESTIMATOR");

  if (model.is_quantal())
    boost::throw_exception(std::invalid_argument("not classical Ising model"));
  if (model.has_field()) std::cerr << "WARNING: model has magnetic field\n";
  if (model.is_frustrated())
    std::cerr << "WARNING: model is classically frustrated\n";
  if (!use_improved_estimator)
    std::cerr << "WARNING: improved estimator is disabled\n";

  lattice.generate_virtual_graph(model, model.has_d_term());
  perm.resize(max_virtual_sites(lattice));

  coupling.resize(num_bonds(lattice.vg()));
  prob.resize(num_bonds(lattice.vg()), 2);
  int b = 0;
  looper::real_bond_iterator<lattice_t>::type rbi, rbi_end;
  for (boost::tie(rbi, rbi_end) = bonds(lattice.rg()); rbi != rbi_end; ++rbi) {
    double jz = -model.bond(*rbi, lattice.rg()).jz;
      // positive for ferromagnetic
    looper::virtual_bond_iterator<lattice_t>::type vbi, vbi_end;
    for (boost::tie(vbi, vbi_end) = bonds(lattice, *rbi); vbi != vbi_end;
         ++vbi, ++b)
      coupling[b] = jz;
  }
  if (model.has_field()) {
    field.resize(num_sites(lattice.vg()));
    int i = 0;
    looper::real_site_iterator<lattice_t>::type rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(lattice.rg()); rsi != rsi_end;
         ++rsi) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) = sites(lattice, *rsi); vsi != vsi_end;
           ++vsi, ++i)
        field[i] = model.site(*rsi, lattice.rg()).hz;
    }
  }

  // configuration
  int nvs = num_sites(lattice.vg());
  spins.resize(nvs); std::fill(spins.begin(), spins.end(), 0 /* all up */);

  // measurements
  obs << make_observable(alps::SimpleRealObservable("Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Inverse Temperature"));
  obs << make_observable(alps::SimpleRealObservable("Number of Sites"));
  obs << make_observable(alps::SimpleRealObservable("Number of Clusters"));
  looper::energy_estimator::initialize(obs, false);
  estimator.initialize(obs, p, lattice, false, use_improved_estimator);
}

void loop_worker::dostep()
{
  if (!mcs.can_work()) return;
  ++mcs;
  beta = 1.0 / temperature(mcs());

  build();

  //   FIELD               IMPROVE
  flip<boost::mpl::true_,  boost::mpl::true_ >();
  flip<boost::mpl::true_,  boost::mpl::false_>();
  flip<boost::mpl::false_, boost::mpl::true_ >();
  flip<boost::mpl::false_, boost::mpl::false_>();

  measure();
}


//
// cluster construction
//

void loop_worker::build()
{
  // initialize cluster information (setup cluster fragments)
  int nvs = num_sites(lattice.vg());
  fragments.resize(0); fragments.resize(nvs);

  // build table of probabilities
  int nvb = num_bonds(lattice.vg());
  if (mcs() <= temperature.annealing_steps()+1)
    for (int b = 0; b < nvb; ++b) {
      prob(b, 0) = 1-std::exp(-0.5*beta*coupling[b]); // for parallel
      prob(b, 1) = 1-std::exp( 0.5*beta*coupling[b]); // for antiparallel
    }

  // building up clusters
  for (int b = 0; b < nvb; ++b) {
    int s0 = source(b, lattice.vg());
    int s1 = target(b, lattice.vg());
    if (random() < prob(b, spins[s0] ^ spins[s1])) unify(fragments, s0, s1);
  }

  // symmetrize spins
  if (max_virtual_sites(lattice) > 1) {
    looper::real_site_iterator<lattice_t>::type rsi, rsi_end;
    for (boost::tie(rsi, rsi_end) = sites(lattice.rg()); rsi != rsi_end;
         ++rsi) {
      looper::virtual_site_iterator<lattice_t>::type vsi, vsi_end;
      boost::tie(vsi, vsi_end) = sites(lattice, *rsi);
      int offset = *vsi;
      int s2 = *vsi_end - *vsi;
      for (int i = 0; i < s2; ++i) perm[i] = i;
      looper::partitioned_random_shuffle(perm.begin(), perm.begin() + s2,
        spins.begin() + offset, spins.begin() + offset, random);
      for (int i = 0; i < s2; ++i)
        unify(fragments, offset+i, offset+perm[i]);
    }
  }
}


//
// cluster flip
//

template<typename FIELD, typename IMPROVE>
void loop_worker::flip()
{
  if (!(false == FIELD() && use_improved_estimator == IMPROVE())) return;

  int nvs = num_sites(lattice.vg());

  // assign cluster id
  int nc = 0;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) if (fi->is_root()) fi->id = nc++;
  for (std::vector<cluster_fragment_t>::iterator fi = fragments.begin();
       fi != fragments.end(); ++fi) fi->id = cluster_id(fragments, *fi);
  clusters.resize(0); clusters.resize(nc);

  cluster_info_t::accumulator<cluster_fragment_t, FIELD,
    boost::mpl::false_, IMPROVE> weight(clusters, fragments, field, 0, 0);
  looper::accumulator<estimator_t, cluster_fragment_t, IMPROVE>
    accum(estimates, nc, lattice, estimator, fragments);
  for (unsigned int s = 0; s < nvs; ++s) {
    weight.at_bot(s, time_t(0), s, spins[s]);
    weight.at_top(s, time_t(1), s, spins[s]);
    accum.at_bot(s, time_t(0), s, spins[s]);
    accum.at_top(s, time_t(1), s, spins[s]);
  }

  // determine whether clusters are flipped or not
  for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
       ci != clusters.end(); ++ci) {
    ci->to_flip = ((2*random()-1) < (FIELD() ? std::tanh(ci->weight) : 0));
  }

  // improved measurement
  if (IMPROVE()) {
    typename looper::collector<estimator_t>::type
      coll = get_collector(estimator);
    coll = std::accumulate(estimates.begin(), estimates.end(), coll);
    estimator.improved_measurement(obs, lattice, beta, 1,
      spins, operator_string_t(), spins, fragments, coll);
  }
  obs["Number of Clusters"] << (double)clusters.size();

  // flip spins
  for (int s = 0; s < nvs; ++s)
    if (clusters[fragments[s].id].to_flip) spins[s] ^= 1;
}


//
// measurement
//

void loop_worker::measure()
{
  int nrs = num_sites(lattice.rg());
  obs["Temperature"] << 1/beta;
  obs["Inverse Temperature"] << beta;
  obs["Number of Sites"] << (double)nrs;

  // energy
  double ene = energy_offset;
  int nvb = num_bonds(lattice.vg());
  for (int b = 0; b < nvb; ++b) {
    int s0 = source(b, lattice.vg());
    int s1 = target(b, lattice.vg());
    ene -= 0.5 * coupling[b] * (0.5 - (spins[s0] ^ spins[s1]));
  }
  if (model.has_field()) {
    for (std::vector<cluster_info_t>::iterator ci = clusters.begin();
         ci != clusters.end(); ++ci)
      if (ci->to_flip) ene -= ci->weight;
      else ene += ci->weight;
  }
  looper::energy_estimator::measurement(obs, lattice, beta, 0, 1, ene);

  estimator.normal_measurement(obs, lattice, beta, 1, spins,
                               operator_string_t(), spins);
}


//
// dynamic registration to the factories
//

const bool loop_registered =
  loop_factory::instance()->register_worker<loop_worker>("Ising");
const bool evaluator_registered = evaluator_factory::instance()->
  register_evaluator<looper::evaluator<loop_config::measurement_set> >("Ising");

} // end namespace
