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

#include "observable.h"
#include "options.h"
#include <looper/union_find.h>
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
  cluster_t(bool t = false) : to_flip(t), size(0), mag(0), length(0) {}
  bool to_flip;
  int size;
  int mag;
  double length;
};

typedef looper::union_find::node fragment_t;

// lattice helper functions (returns site index at left/right end of a bond)
inline int left(int /* L */, int b) { return b; }
inline int right(int L, int b) { return (b == L-1) ? 0 : b+1; }

// helper function, which returns x^2
template<typename T> inline T power2(T x) { return x * x; }

int main(int argc, char* argv[]) {
  // parameters
  options p(argc, argv, true, false);
  const unsigned int nsites = p.length;
  const unsigned int nbonds = nsites;
  const unsigned int sweeps = p.sweeps;
  const unsigned int therm = p.therm;
  const double beta = 1. / p.temperature;

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
  std::vector<int> spins_c(nsites);
  std::vector<int> current(nsites);

  // cluster information
  std::vector<fragment_t> fragments;
  std::vector<cluster_t> clusters;

  // oservables
  observable energy;
  observable smag; // staggered magnetizetion^2
  observable ssus; // staggered susceptibility
  observable usus; // uniform susceptibility

  //
  // Monte Carlo steps
  //

  boost::timer tm;

  for (unsigned int mcs = 0; mcs < therm + sweeps; ++mcs) {

    //
    // diagonal update and cluster construction
    //

    // initialize spin & operator information
    std::copy(spins.begin(), spins.end(), spins_c.begin());
    std::swap(operators, operators_p); operators.resize(0);

    // initialize cluster information (setup s cluster fragments)
    fragments.resize(0); fragments.resize(nsites);
    for (int s = 0; s < nsites; ++s) current[s] = s;

    double t = r_time();
    for (std::vector<local_operator_t>::iterator opi = operators_p.begin();
         t < 1 || opi != operators_p.end();) {

      // diagonal update
      if (opi == operators_p.end() || t < opi->time) {
        int b = nbonds * r_uniform();
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

    for (int s = 0; s < nsites; ++s) unify(fragments, s, current[s]);

    //
    // cluster flip
    //

    // assign cluster id
    int nc = 0;
    BOOST_FOREACH(fragment_t& f, fragments) if (f.is_root()) f.set_id(nc++);
    BOOST_FOREACH(fragment_t& f, fragments) f.set_id(cluster_id(fragments, f));
    clusters.resize(0); clusters.resize(nc);

    // 'flip' operators & do improved measurements
    BOOST_FOREACH(local_operator_t& op, operators) {
      double t = op.time;
      clusters[fragments[op.lower_loop].id()].length += 2 * t;
      clusters[fragments[op.upper_loop].id()].length -= 2 * t;
    }
    for (unsigned int s = 0; s < nsites; ++s) {
      int id = fragments[s].id();
      clusters[id].size += 1;
      clusters[id].mag += 1 - 2 * spins[s];
      clusters[id].length += 1;
    }

    // determine whether clusters are flipped or not
    BOOST_FOREACH(cluster_t& ci, clusters) ci.to_flip = (r_uniform() < 0.5);

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

    // accumurate loop length and magnetization
    double s2 = 0;
    double m2 = 0;
    double l2 = 0;
    for (std::vector<cluster_t>::const_iterator pi = clusters.begin();
         pi != clusters.end(); ++pi) {
      s2 += power2(pi->size);
      m2 += power2(pi->mag);
      l2 += power2(pi->length);
    }

    energy << (0.25 * nbonds - operators.size() / beta) / nsites;
    smag << 0.25 * s2 / nsites;
    usus << 0.25 * beta * m2 / nsites;
    ssus << 0.25 * beta * l2 / nsites;
  }

  std::cerr << "Speed = " << (therm + sweeps) / tm.elapsed()
            << " MCS/sec\n";
  std::cout << "Energy per Site           = "
            << energy.mean() << " +- " << energy.error() << std::endl
            << "Staggered Magnetization^2 = "
            << smag.mean() << " +- " << smag.error() << std::endl
            << "Uniform Susceptibility    = "
            << usus.mean() << " +- " << usus.error() << std::endl
            << "Staggered Susceptibility  = "
            << ssus.mean() << " +- " << ssus.error() << std::endl;
}
