/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

// Spin-1/2 Antiferromagnetic Heisenberg Chain
// [continuous time path integral; using std::vector<> for operator string]

#include "common.h"
#include "observable.h"
#include "options.h"

#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>

#include <algorithm> // for std::swap
#include <cstdlib> // for std::exit
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
  // parameters
  options p(argc, argv);
  if (!p.valid) std::exit(-1);
  const unsigned int nsites = p.length;
  const unsigned int nbonds = nsites;
  const unsigned int sweeps = p.sweeps;
  const unsigned int therm = p.therm;
  const double beta = 1 / p.temperature;

  // random number generators
  boost::mt19937 eng(29833u);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(beta * nbonds / 2));

  // vector of operators
  std::vector<local_operator_t> operators, operators_p;

  // spin configuration at t = 0 (1 for down and 0 for up)
  std::vector<int> spins(nsites);
  std::fill(spins.begin(), spins.end(), 0 /* all up */);
  std::vector<int> current(nsites);

  // cluster information
  std::vector<fragment_t> fragments;
  std::vector<bool> to_flip;
  std::vector<estimate_t> estimates;

  // oservables
  observable num_clusters;
  observable energy;
  observable usus; // uniform susceptibility
  observable smag; // staggered magnetizetion^2
  observable ssus; // staggered susceptibility

  //
  // Monte Carlo steps
  //

  boost::timer tm;

  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {

    //
    // diagonal update and cluster construction
    //

    // initialize operator information
    std::swap(operators, operators_p); operators.resize(0);

    // initialize cluster information
    fragments.resize(0); fragments.resize(nsites);
    for (int s = 0; s < nsites; ++s) current[s] = s;

    double t = r_time();
    for (std::vector<local_operator_t>::iterator opi = operators_p.begin();
         t < 1 || opi != operators_p.end();) {

      // diagonal update
      if (opi == operators_p.end() || t < opi->time) {
        const int b = static_cast<int>(nbonds * r_uniform());
        if (spins[left(nbonds, b)] != spins[right(nbonds, b)]) {
          operators.push_back(local_operator_t(b, t));
          t += r_time();
        } else {
          t += r_time();
          continue;
        }
      } else {
        if (opi->type == diagonal) {
          ++opi;
          continue;
        } else {
          operators.push_back(*opi);
          ++opi;
        }
      }

      std::vector<local_operator_t>::iterator oi = operators.end() - 1;
      const int s0 = left(nbonds, oi->bond);
      const int s1 = right(nbonds, oi->bond);
      if (oi->type == offdiagonal) {
        spins[s0] ^= 1;
        spins[s1] ^= 1;
      }
      oi->lower_cluster = unify(fragments, current[s0], current[s1]);
      oi->upper_cluster = current[s0] = current[s1] = add(fragments);
    }

    for (int s = 0; s < nsites; ++s) unify(fragments, s, current[s]);

    //
    // cluster flip
    //

    // assign cluster id
    int nc = 0;
    BOOST_FOREACH(fragment_t& f, fragments) if (f.is_root()) f.set_id(nc++);
    BOOST_FOREACH(fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
    to_flip.resize(nc);
    estimates.resize(0); estimates.resize(nc);

    BOOST_FOREACH(local_operator_t& op, operators) {
      const double t = op.time;
      estimates[fragments[op.lower_cluster].id()].length += 2 * t;
      estimates[fragments[op.upper_cluster].id()].length -= 2 * t;
    }
    for (unsigned int s = 0; s < nsites; ++s) {
      const int id = fragments[s].id();
      estimates[id].mag += 1 - 2 * spins[s];
      estimates[id].size += 1;
      estimates[id].length += 1;
    }

    // accumulate cluster properties
    collector_t coll;
    coll.set_num_operators(operators.size());
    coll.set_num_clusters(nc);
    BOOST_FOREACH(estimate_t const& est, estimates) coll += est;

    // determine whether clusters are flipped or not
    for (int c = 0; c < nc; ++c) to_flip[c] = (r_uniform() < 0.5);

    // flip operators & spins
    BOOST_FOREACH(local_operator_t& op, operators)
      if (to_flip[fragments[op.lower_cluster].id()] ^
          to_flip[fragments[op.upper_cluster].id()]) op.flip();
    for (int s = 0; s < nsites; ++s)
      if (to_flip[fragments[s].id()]) spins[s] ^= 1;

    //
    // measurements
    //

    if (mcs >= therm) {
      num_clusters << coll.num_clusters();
      energy << (0.25 * nbonds - coll.num_operators() / beta) / nsites;
      usus << 0.25 * beta * coll.usus / nsites;
      smag << 0.25 * coll.smag;
      ssus << 0.25 * beta * coll.ssus / nsites;
    }
  }

  std::cerr << "Speed = " << (therm + sweeps) / tm.elapsed()
            << " MCS/sec\n";
  std::cout << "Number of Clusters        = "
            << num_clusters.mean() << " +- " << num_clusters.error() << std::endl
            << "Energy Density            = "
            << energy.mean() << " +- " << energy.error() << std::endl
            << "Uniform Susceptibility    = "
            << usus.mean() << " +- " << usus.error() << std::endl
            << "Staggered Magnetization^2 = "
            << smag.mean() << " +- " << smag.error() << std::endl
            << "Staggered Susceptibility  = "
            << ssus.mean() << " +- " << ssus.error() << std::endl;
}
