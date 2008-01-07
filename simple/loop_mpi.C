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

// Spin-1/2 Antiferromagnetic Heisenberg Chain (parallel version)
// [continuous time path integral; using std::vector<> for operator string]

#include "observable.h"
#include "options.h"
#include <looper/union_find.h>
#include <looper/parallel.h>
#include <algorithm> // for std::swap
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <iostream>
#include <vector>

enum operator_type { diagonal, offdiagonal };

struct local_operator_t {
  local_operator_t() {}
  local_operator_t(int b, double t) : type(diagonal), bond(b), time(t) {}
  void flip() { type = (type == diagonal ? offdiagonal : diagonal); }
  operator_type type;
  unsigned int bond;
  unsigned int upper_loop, lower_loop;
  double time;
};

struct cluster_t {
  cluster_t(bool t = false) : to_flip(t) {}
  bool to_flip;
};

struct estimate_t {
  estimate_t() : mag(0), size(0), length(0) {}
  double mag;
  double size;
  double length;
};

struct accumulate_t {
  accumulate_t() : nop(0), usus(0), smag(0), ssus(0) {}
  void add(accumulate_t const& accum) {
    nop += accum.nop;
    usus += accum.usus;
    smag += accum.smag;
    ssus += accum.ssus;
  }
  void add_cluster(estimate_t const& est) {
    usus += est.mag * est.mag;
    smag += est.size * est.size;
    ssus += est.length * est.length;
  }
  double nop;
  double usus;
  double smag;
  double ssus;
};

typedef looper::union_find::node fragment_t;

// lattice helper functions (returns site index at left/right end of a bond)
inline int left(int /* L */, int b) { return b; }
inline int right(int L, int b) { return (b == L-1) ? 0 : b+1; }

// helper function, which returns x^2
template<typename T> inline T power2(T x) { return x * x; }

int main(int argc, char* argv[]) {

  MPI_Init(&argc, &argv);
  int num_processes, process_id;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

  // parameters
  options p(argc, argv, true, false, process_id == 0);
  const unsigned int nsites = p.length;
  const unsigned int nbonds = nsites;
  const unsigned int sweeps = p.sweeps;
  const unsigned int therm = p.therm;
  const double beta = 1. / p.temperature;
  const double tau0 = 1. * process_id / num_processes;
  const double tau1 = 1. * (process_id+1) / num_processes;

  // random number generators
  boost::mt19937 eng(29833u & (process_id << 8));
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    r_uniform(eng, boost::uniform_real<>());
  boost::variate_generator<boost::mt19937&, boost::exponential_distribution<> >
    r_time(eng, boost::exponential_distribution<>(beta * nbonds / 2));

  // vector of operators
  std::vector<local_operator_t> operators, operators_p;

  // spin configuration at t = tau0 (1 for down and 0 for up)
  std::vector<int> spins(nsites);
  std::fill(spins.begin(), spins.end(), 0 /* all up */);
  std::vector<int> spins_c(nsites);
  std::vector<int> current(nsites);

  // cluster information
  std::vector<fragment_t> fragments;
  std::vector<cluster_t> clusters;
  std::vector<estimate_t> estimates;

  // oservables
  observable energy;
  observable usus; // uniform susceptibility
  observable smag; // staggered magnetizetion^2
  observable ssus; // staggered susceptibility

  // helper for parallelization
  typedef looper::parallel_cluster_unifier<fragment_t, estimate_t, accumulate_t> unifier_t;
  unifier_t unifier(MPI_COMM_WORLD, nsites);

  //
  // Monte Carlo steps
  //

  MPI_Barrier(MPI_COMM_WORLD);
  boost::timer tm;

  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {

    //
    // diagonal update and cluster construction
    //

    // initialize spin & operator information
    std::copy(spins.begin(), spins.end(), spins_c.begin());
    std::swap(operators, operators_p); operators.resize(0);

    // initialize cluster information (setup s cluster fragments)
    fragments.resize(0); fragments.resize(2 * nsites);
    for (int s = 0; s < nsites; ++s) current[s] = s;

    double t = tau0 + r_time();
    for (std::vector<local_operator_t>::iterator opi = operators_p.begin();
         t < tau1 || opi != operators_p.end();) {

      // diagonal update
      if (opi == operators_p.end() || t < opi->time) {
        int b = static_cast<int>(nbonds * r_uniform());
        if (spins_c[left(nbonds, b)] != spins_c[right(nbonds, b)]) {
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
      int s0 = left(nbonds, oi->bond);
      int s1 = right(nbonds, oi->bond);
      if (oi->type == offdiagonal) {
        spins_c[s0] ^= 1;
        spins_c[s1] ^= 1;
      }
      oi->lower_loop = unify(fragments, current[s0], current[s1]);
      oi->upper_loop = current[s0] = current[s1] = add(fragments);
    }

    for (int s = 0; s < nsites; ++s) unify(fragments, current[s], nsites + s);
    for (int s = 0; s < 2 * nsites; ++s) {
      int r = root_index(fragments, s);
      if (r > s) {
        fragments[s].set_as_root(fragments[r].weight());
        fragments[r].set_parent(s);
      }
    }

    //
    // cluster flip
    //

    // assign cluster id
    int nc = 0;
    for (int s = 0; s < 2 * nsites; ++s) if (fragments[s].is_root()) fragments[s].set_id(nc++);
    int noc = nc;
    for (int s = 2 * nsites; s < fragments.size(); ++s)
      if (fragments[s].is_root()) fragments[s].set_id(nc++);
    BOOST_FOREACH(fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
    clusters.resize(0); clusters.resize(nc);
    estimates.resize(0); estimates.resize(nc);

    for (unsigned int s = 0; s < nsites; ++s) {
      int id = fragments[s].id();
      estimates[id].mag += 1 - 2 * spins[s];
      estimates[id].size -= tau0;
      estimates[id].length -= tau0;
    }
    BOOST_FOREACH(local_operator_t& op, operators) {
      double t = op.time;
      estimates[fragments[op.lower_loop].id()].length += 2 * t;
      estimates[fragments[op.upper_loop].id()].length -= 2 * t;
    }
    for (unsigned int s = 0; s < nsites; ++s) {
      int id = fragments[s].id();
      estimates[id].size += tau1;
      estimates[id].length += tau1;
    }

    // accumulate loop length and magnetization
    accumulate_t accum;
    accum.nop = operators.size();
    for (int c = noc; c < nc; ++c) accum.add_cluster(estimates[c]);

    // global unification of open clusters
    std::vector<unifier_t::flip_t> const& to_flip =
      unifier.unify(noc, fragments, estimates, accum, r_uniform);

    // determine whether clusters are flipped or not
    for (int c = 0; c < noc; ++c) clusters[c].to_flip = to_flip[c].flip();
    for (int c = noc; c < nc; ++c) clusters[c].to_flip = (r_uniform() < 0.5);

    // flip operators & spins
    BOOST_FOREACH(local_operator_t& op, operators)
      if (clusters[fragments[op.lower_loop].id()].to_flip ^
          clusters[fragments[op.upper_loop].id()].to_flip) op.flip();
    for (int s = 0; s < nsites; ++s)
      if (clusters[fragments[s].id()].to_flip) spins[s] ^= 1;

    if (mcs < therm) continue;

    //
    // measurements
    //

    if (process_id == 0) {
      energy << (0.25 * nbonds - accum.nop / beta) / nsites;
      usus << 0.25 * beta * accum.usus / nsites;
      smag << 0.25 * accum.smag / nsites;
      ssus << 0.25 * beta * accum.ssus / nsites;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (process_id == 0) {
    std::cerr << "Speed = " << (therm + sweeps) / tm.elapsed()
              << " MCS/sec\n";
    std::cout << "Energy per Site           = "
              << energy.mean() << " +- " << energy.error() << std::endl
              << "Uniform Susceptibility    = "
              << usus.mean() << " +- " << usus.error() << std::endl
              << "Staggered Magnetization^2 = "
              << smag.mean() << " +- " << smag.error() << std::endl
              << "Staggered Susceptibility  = "
              << ssus.mean() << " +- " << ssus.error() << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
