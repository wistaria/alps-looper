/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>
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

#include "union_find.h"

enum operator_type { diagonal, offdiagonal };

struct local_operator_t {
  local_operator_t() {}
  local_operator_t(unsigned int b, double t) : type(diagonal), bond(b), time(t) {}
  void flip() { type = (type == diagonal ? offdiagonal : diagonal); }
  operator_type type;
  unsigned int bond;
  unsigned int upper_cluster, lower_cluster;
  double time;
};

struct cluster_t {
  cluster_t(bool t = false) : to_flip(t) {}
  bool to_flip;
};

struct estimate_t {
  estimate_t() : mag(0), size(0), length(0) {}
  estimate_t& operator+=(estimate_t const& est) {
    mag += est.mag;
    size += est.size;
    length += est.length;
    return *this;
  }
  double mag;
  double size;
  double length;
};

struct collector_t {
  collector_t() : noc(0), nc(0), nop(0), usus(0), smag(0), ssus(0) {}
  collector_t& operator+=(collector_t const& coll) {
    nc += coll.nc;
    nop += coll.nop;
    usus += coll.usus;
    smag += coll.smag;
    ssus += coll.ssus;
    return *this;
  }
  collector_t& operator+= (estimate_t const& est) {
    usus += est.mag * est.mag;
    smag += est.size * est.size;
    ssus += est.length * est.length;
    return *this;
  }
  void set_num_open_clusters(unsigned int n) { noc += n; }
  unsigned int num_open_clusters() const { return noc; }
  void set_num_clusters(unsigned int n) { nc = n; }
  void inc_num_clusters(unsigned int n) { nc += n; }
  double num_clusters() const { return nc; }
  void set_num_operators(unsigned int n) { nop = n; }
  double num_operators() const { return nop; }

  unsigned int noc; // number of open clusters (for parallel QMC)
  double nc;        // total number of (closed) clusters
  double nop;       // total number of operators
  double usus;
  double smag;
  double ssus;
};

typedef looper::union_find::node fragment_t;

// lattice helper functions (returns site index at left/right end of a bond)
inline int left(unsigned int /* L */, unsigned int b) { return b; }
inline int right(unsigned int L, unsigned int b) { return (b == L-1) ? 0 : b+1; }
