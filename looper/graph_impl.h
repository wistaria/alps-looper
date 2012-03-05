/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_GRAPH_IMPL_H
#define LOOPER_GRAPH_IMPL_H

#include <alps/config.h>
#if defined(LOOPER_ENABLE_OPENMP) && defined(ALPS_ENABLE_OPENMP_WORKER) && !defined(LOOPER_OPENMP)
# define LOOPER_OPENMP
#endif

#include "graph.h"
#include "location_impl.h"
#include "random_choice.h"
#ifdef HAVE_PARAPACK_13
# include <alps/math.hpp>
#else
# include <alps/numeric/is_nonzero.hpp>
# include <alps/numeric/is_zero.hpp>
#endif
#include <boost/array.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <stdexcept>

#ifdef LOOPER_OPENMP
# include <omp.h>
#endif

#ifdef HAVE_PARAPACK_13
using alps::is_nonzero;
using alps::is_zero;
#else
using alps::numeric::is_nonzero;
using alps::numeric::is_zero;
#endif

namespace looper {

//
// site graph type
//
//   g = 0 (vertically disconnected)
//

struct site_graph_type {
  static bool is_valid_gid(int g) { return g == 0; }
  static bool is_compatible(int /* g */, int /* c */) { return true; }
  static bool is_frozen(int /* g */) { return false; }
  template<class T>
  // obsolete interface
  static boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(int /* g */, std::vector<T>& fragments, int curr) {
    int loop1 = add(fragments);
    return boost::make_tuple(loop1, curr, loop1);
  }
  // current interface
  template<class T>
  static boost::tuple<int /* fid */, int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(int /* g */, std::vector<T>& fragments, int fid, int curr) {
    fragments[fid] = T();
    int loop1 = fid;
    ++fid;
    return boost::make_tuple(fid, loop1, curr, loop1);
  }
};


//
// bond graph types
//
//   g = 0   o.o    g = 1   o o
//                           X
//           o.o            o o
//
//   g = 2   o.o    g = 3   o-o
//           |*|            |X|
//           o.o            o-o
//
//   g = 4   o o    g = 5   o-o
//            *
//           o o            o-o

// optimized bond_graph_type for Ising model

struct ising_bond_graph_type {
  static bool is_valid_gid(int g) { return g == 2 || g == 3; }
  static bool is_compatible(int g, int c0, int c1) { return (g & 1) ^ c0 ^ c1; }
  static bool is_frozen(int /* g */) { return true; }
  // obsolete interface
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop = unify(fragments, curr0, curr1);
    return boost::make_tuple(loop, loop, loop, loop);
  }
  // current interface
  template<class T>
  static boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */,
                      int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int fid, int curr0, int curr1) {
    int loop = unify(fragments, curr0, curr1);
    return boost::make_tuple(fid, loop, loop, loop, loop);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return false; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_diagonal(RNG& rng, std::vector<diagonal_choice_helper> const&, int,
    int c0, int c1) {
    return (c0 == c1) ? 3 : 2;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return false; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG&, std::vector<offdiagonal_choice_helper> const&, int, int,
    int) {
    boost::throw_exception(std::logic_error("ising_bond_graph_type::choose_offdiagonal"));
    return 0;
  }
};

// optimized bond_graph_type for Heisenberg Antiferromagnet

struct haf_bond_graph_type {
  static bool is_valid_gid(int g) { return g == 0; }
  static bool is_compatible(int /* g */, int c0, int c1) { return c0 != c1; }
  static bool is_frozen(int /* g */) { return false; }
  // obsolete interface
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0 = unify(fragments, curr0, curr1);
    int loop1 = curr0 = curr1 = add(fragments);
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }
  // current interface
  template<class T>
  static boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */,
                      int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int fid, int curr0, int curr1) {
    fragments[fid] = T();
    int loop0 = unify(fragments, curr0, curr1);
    int loop1 = curr0 = curr1 = fid;
    ++fid;
    return boost::make_tuple(fid, curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return false; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_diagonal(RNG&, std::vector<diagonal_choice_helper> const&, int, int, int) {
    return 0;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return false; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG&, std::vector<offdiagonal_choice_helper> const&, int, int,
    int) {
    return 0;
  }
};

// optimized bond_graph_type for Heisenberg Ferromagnet

struct hf_bond_graph_type {
  static bool is_valid_gid(int g) { return g == 1; }
  static bool is_compatible(int c0, int c1) { return c0 == c1; }
  static bool is_frozen(int /* g */) { return false; }
  // obsolete interface
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0 = curr0;
    int loop1 = curr1;
    std::swap(curr0, curr1);
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }
  // current interface
  template<class T>
  static boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */,
                      int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int fid, int curr0, int curr1) {
    int loop0 = curr0;
    int loop1 = curr1;
    std::swap(curr0, curr1);
    return boost::make_tuple(fid, curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return false; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_diagonal(RNG&, std::vector<diagonal_choice_helper> const&, int, int, int) {
    return 1;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return false; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const&) {}
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG&, std::vector<offdiagonal_choice_helper> const&, int, int,
    int) {
    return 1;
  }
};

// bond_graph_type for XXZ interaction

struct xxz_bond_graph_type {
  static bool is_valid_gid(int g) { return g >= 0 && g <= 3; }
  static bool is_compatible(int g, int c0, int c1) { return (g & 1) ^ c0 ^ c1; }
  static bool is_frozen(int g) { return ((g & 2) == 2); }
  // obsolete interface
  template<class T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0, loop1;
    if ((g & 2) == 2) {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    } else if (g == 0) {
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = add(fragments);
    } else {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    }
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }
  // current interface
  template<class T>
  static boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */,
                      int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int fid, int curr0, int curr1) {
    int loop0, loop1;
    if ((g & 2) == 2) {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    } else if (g == 0) {
      fragments[fid] = T();
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = fid;
      ++fid;
    } else {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    }
    return boost::make_tuple(fid, curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return true; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const& v) {
      p[0] = is_nonzero<1>(v[0] + v[2]) ? v[0] / (v[0] + v[2]) : 1; // for antiparallel
      p[1] = is_nonzero<1>(v[1] + v[3]) ? v[1] / (v[1] + v[3]) : 1; // for parallel
    }
    boost::array<double, 2> p;
  };
  template<typename RNG>
  static int choose_diagonal(RNG& rng, std::vector<diagonal_choice_helper> const& helpers, int b,
    int c0, int c1) {
    int c = 1 ^ c0 ^ c1; // 0 for antiparallel, 1 for parallel,
    return (rng() < helpers[b].p[c] ? 0 : 2) ^ c;
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return true; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const& v) {
      p = is_nonzero<1>(v[0] + v[1]) ? v[0] / (v[0] + v[1]) : 1;
    }
    double p;
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG& rng, std::vector<offdiagonal_choice_helper> const& helpers,
    int b, int /* c0 */, int /* c1 */) {
    return rng() < helpers[b].p ? 0 : 1;
  }
};

// bond_graph_type for XYZ interaction

struct xyz_bond_graph_type {
  static bool is_valid_gid(int g) { return g >= 0 && g <= 5; }
  static bool is_compatible(int g, int c0, int c1) { return (g & 1) ^ c0 ^ c1; }
  static bool is_frozen(int g) { return ((g & 2) == 2); }
  // obsolete interface
  template<typename T>
  static boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int curr0, int curr1) {
    int loop0, loop1;
    if ((g & 2) == 2) {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    } else if (g == 0 || g == 5) {
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = add(fragments);
    } else {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    }
    return boost::make_tuple(curr0, curr1, loop0, loop1);
  }
  // current interface
  template<typename T>
  static boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */,
                      int /* loop1 */>
  reconnect(int g, std::vector<T>& fragments, int fid, int curr0, int curr1) {
    int loop0, loop1;
    if ((g & 2) == 2) {
      loop0 = loop1 = curr0 = curr1 = unify(fragments, curr0, curr1);
    } else if (g == 0 || g == 5) {
      fragments[fid] = T();
      loop0 = unify(fragments, curr0, curr1);
      loop1 = curr0 = curr1 = fid;
      ++fid;
    } else {
      loop0 = curr0;
      loop1 = curr1;
      std::swap(curr0, curr1);
    }
    return boost::make_tuple(fid, curr0, curr1, loop0, loop1);
  }

  static bool has_nontrivial_diagonal_choice_helper() { return true; }
  struct diagonal_choice_helper {
    diagonal_choice_helper() {}
    template<std::size_t N>
    diagonal_choice_helper(boost::array<double, N> const& v) {
      // for antiparallel
      p[0] = is_nonzero<1>(v[0] + v[2] + v[4]) ? v[0] / (v[0] + v[2] + v[4]) : 1;
      p[2] = is_nonzero<1>(v[0] + v[2] + v[4]) ? (v[0] + v[2])/ (v[0] + v[2] + v[4]) : 1;
      // for parallel
      p[1] = is_nonzero<1>(v[1] + v[3] + v[5]) ? v[1] / (v[1] + v[3] + v[5]) : 1;
      p[3] = is_nonzero<1>(v[1] + v[3] + v[5]) ? (v[1] + v[3]) / (v[1] + v[3] + v[5]) : 1;
    }
    boost::array<double, 4> p;
  };
  template<typename RNG>
  static int choose_diagonal(RNG& rng, std::vector<diagonal_choice_helper> const& helpers, int b,
    int c0, int c1) {
    double r = rng();
    if (c0 != c1) {
      if (r < helpers[b].p[0])
        return 0;
      else if (r < helpers[b].p[2])
        return 2;
      else
        return 4;
    } else {
      if (r < helpers[b].p[1])
        return 1;
      else if (r < helpers[b].p[3])
        return 3;
      else
        return 5;
    }
  }

  static bool has_nontrivial_offdiagonal_choice_helper() { return true; }
  struct offdiagonal_choice_helper {
    offdiagonal_choice_helper() {}
    template<std::size_t N>
    offdiagonal_choice_helper(boost::array<double, N> const& v) {
      p[0] = is_nonzero<1>(v[0] + v[1]) ? v[0] / (v[0] + v[1]) : 1;
      p[1] = is_nonzero<1>(v[4] + v[5]) ? v[4] / (v[4] + v[5]) : 1;
    }
    boost::array<double, 2> p;
  };
  template<typename RNG>
  static int choose_offdiagonal(RNG& rng, std::vector<offdiagonal_choice_helper> const& helpers,
    int b, int c0, int c1) {
    double r = rng();
    if (c0 != c1)
      return (r < helpers[b].p[0] ? 0 : 1);
    else
      return (r < helpers[b].p[1] ? 4 : 5);
  }
};

template<typename SITE, typename BOND, typename LOC>
class local_graph {
public:
  typedef SITE site_graph_t;
  typedef BOND bond_graph_t;
  typedef LOC location_t;

  local_graph() {}
  local_graph(int type, const location_t& loc) : type_(type), loc_(loc) {}

  int type() const { return type_; }
  bool is_compatible(int c) const {
#ifndef NDEBUG
    if (!is_site())
      boost::throw_exception(std::logic_error("local_graph<>::is_compatible"));
#endif
    return site_graph_t::is_compatible(type_, c);
  }
  bool is_compatible(int c0, int c1) const {
#ifndef NDEBUG
    if (!is_bond())
      boost::throw_exception(std::logic_error("local_graph<>::is_compatible"));
#endif
    return bond_graph_t::is_compatible(type_, c0, c1);
  }

  const location_t& loc() const { return loc_; }
  int pos() const { return loc_.pos(); }
  int source() const { return loc_.source(); }
  int target() const { return loc_.target(); }
  bool is_bond() const { return loc_.is_bond(); }
  bool is_site() const { return loc_.is_site(); }

  // obsolete interface
  template<class T>
  boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr) const {
    return site_graph_t::reconnect(type_, fragments, curr);
  }
  // current interface
  template<class T>
  boost::tuple<int /* fid */, int /* curr */, int /* loop0 */, int /* loop1 */>
  reconnect_site(std::vector<T>& fragments, int fid, int curr) const {
    return site_graph_t::reconnect(type_, fragments, fid, curr);
  }

  // obsolete interface
  template<class T>
  boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect(std::vector<T>& fragments, int curr0, int curr1) const {
    return bond_graph_t::reconnect(type_, fragments, curr0, curr1);
  }
  // current interface
  template<class T>
  boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
  reconnect_bond(std::vector<T>& fragments, int fid, int curr0, int curr1) const {
    return bond_graph_t::reconnect(type_, fragments, fid, curr0, curr1);
  }

  static int loop_l(int /* t */, int loop0, int /* loop1 */) { return loop0; }
  static int loop_u(int /* t */, int /* loop0 */, int loop1) { return loop1; }
  static int loop_l0(int /* t */, int loop0, int /* loop1 */) { return loop0; }
  static int loop_l1(int t, int loop0, int loop1) { return (t ==0 || t == 5) ? loop0 : loop1; }
  static int loop_u0(int /* t */, int /* loop0 */, int loop1) { return loop1; }
  static int loop_u1(int t, int loop0, int loop1) { return (t ==0 || t == 5) ? loop1 : loop0; }

private:
  int type_;
  location_t loc_;
};

template<typename SITE, typename BOND, typename LOC>
inline int type(const local_graph<SITE, BOND, LOC>& g) { return g.type(); }
template<typename SITE, typename BOND, typename LOC>
inline bool is_compatible(const local_graph<SITE, BOND, LOC>& g, int c) {
  return g.is_compatible(c);
}
template<typename SITE, typename BOND, typename LOC>
inline bool is_compatible(const local_graph<SITE, BOND, LOC>& g, int c0, int c1) {
  return g.is_compatible(c0, c1);
}

template<typename SITE, typename BOND, typename LOC>
inline int pos(const local_graph<SITE, BOND, LOC>& g) { return g.pos(); }
template<typename SITE, typename BOND, typename LOC>
inline bool is_site(const local_graph<SITE, BOND, LOC>& g) { return g.is_site(); }
template<typename SITE, typename BOND, typename LOC>
inline bool is_bond(const local_graph<SITE, BOND, LOC>& g) { return g.is_bond(); }

// obsolete interface
template<typename T, typename SITE, typename BOND, typename LOC>
boost::tuple<int /* curr */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph<SITE, BOND, LOC>& g,
          int curr) {
  return g.reconnect(fragments, curr);
}
// current interface
template<typename T, typename SITE, typename BOND, typename LOC>
boost::tuple<int /* fid */, int /* curr */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, int fid, const local_graph<SITE, BOND, LOC>& g,
          int curr) {
  return g.reconnect_site(fragments, fid, curr);
}

// obsolete interface
template<typename T, typename SITE, typename BOND, typename LOC>
boost::tuple<int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, const local_graph<SITE, BOND, LOC>& g,
          int curr0, int curr1) {
  return g.reconnect(fragments, curr0, curr1);
}
// current interface
template<typename T, typename SITE, typename BOND, typename LOC>
boost::tuple<int /* fid */, int /* curr0 */, int /* curr1 */, int /* loop0 */, int /* loop1 */>
reconnect(std::vector<T>& fragments, int fid, const local_graph<SITE, BOND, LOC>& g,
          int curr0, int curr1) {
  return g.reconnect_bond(fragments, fid, curr0, curr1);
}


//
// graph_chooser
//

template<typename LOCAL_GRAPH> class graph_chooser;

template<typename SITE, typename BOND, typename LOC>
class graph_chooser<local_graph<SITE, BOND, LOC> > {
public:
  typedef SITE site_graph_t;
  typedef BOND bond_graph_t;
  typedef LOC location_t;
  typedef local_graph<site_graph_t, bond_graph_t, location_t> local_graph_t;

  graph_chooser() {}
  template<typename WEIGHT_TABLE>
  graph_chooser(const WEIGHT_TABLE& wt, bool is_path_integral) {
    init(wt, is_path_integral);
  }
  template<class WEIGHT_TABLE, typename SHARING>
  graph_chooser(const WEIGHT_TABLE& wt, SHARING const& sharing, bool is_path_integral) {
    init(wt, sharing, is_path_integral);
  }

  template<typename WEIGHT_TABLE>
  void init(const WEIGHT_TABLE& wt, bool is_path_integral) {
    int num_threads = 1;
    dist_graph_.resize(num_threads);
    graph_.resize(num_threads);
    for (int i = 0; i < graph_.size(); ++i) graph_[i].resize(0);
    diag_.resize(0);
    offdiag_.resize(0);
    std::vector<double> w;
    weight_.resize(num_threads);
    for (int i = 0; i < graph_.size(); ++i) weight_[i] = 0;
    BOOST_FOREACH(typename WEIGHT_TABLE::site_weight_t const& sw, wt.site_weights()) {
      for (int g = 0; g < WEIGHT_TABLE::num_site_graphs; ++g) {
        if (is_nonzero<1>(sw.second.v[g])) {
          if (site_graph_t::is_valid_gid(g)) {
            graph_[0].push_back(local_graph_t(g, location_t::site_location(sw.first)));
            w.push_back(sw.second.v[g]);
            weight_[0] += sw.second.v[g];
          } else {
            boost::throw_exception(std::runtime_error("graph_chooser::init"));
          }
        }
      }
    }
    BOOST_FOREACH(typename WEIGHT_TABLE::bond_weight_t const& bw, wt.bond_weights()) {
      for (int g = 0; g < WEIGHT_TABLE::num_bond_graphs; ++g) {
        if (is_nonzero<1>(bw.second.v[g])) {
          if (bond_graph_t::is_valid_gid(g)) {
            graph_[0].push_back(local_graph_t(g, location_t::bond_location(bw.first)));
            w.push_back(bw.second.v[g]);
            weight_[0] += bw.second.v[g];
          } else {
            boost::throw_exception(std::runtime_error("graph_chooser::init"));
          }
        }
      }
      if (!is_path_integral && bond_graph_t::has_nontrivial_diagonal_choice_helper())
        diag_.push_back(typename bond_graph_t::diagonal_choice_helper(bw.second.v));
      if (bond_graph_t::has_nontrivial_offdiagonal_choice_helper())
        offdiag_.push_back(typename bond_graph_t::offdiagonal_choice_helper(bw.second.v));
    }
    if (w.size()) dist_graph_[0].init(w);
    if (is_zero<2>(weight_[0])) weight_[0] = 1.0e-20;
  }

#ifdef LOOPER_OPENMP
  template<typename WEIGHT_TABLE, typename SHARING>
  void init(const WEIGHT_TABLE& wt, SHARING const& sharing, bool is_path_integral) {
    int num_threads = omp_get_max_threads();
    if (num_threads == 1) {
      init(wt, is_path_integral);
      return;
    }
    dist_graph_.resize(num_threads);
    graph_.resize(num_threads);
    diag_.resize(0);
    offdiag_.resize(0);
    weight_.resize(num_threads);
    std::vector<std::vector<double> > w(num_threads);
    // check for site operator
    BOOST_FOREACH(typename WEIGHT_TABLE::site_weight_t const& sw, wt.site_weights()) {
      for (int g = 0; g < WEIGHT_TABLE::num_site_graphs; ++g) {
        if (is_nonzero<1>(sw.second.v[g]))
          boost::throw_exception(std::runtime_error("graph_chooser::init() site operator is not supported in OpenMP parallelization"));
      }
    }
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      graph_[tid].resize(0);
      weight_[tid] = 0;
      // site weight is not supported
      //// BOOST_FOREACH(typename WEIGHT_TABLE::bond_weight_t const& bw, wt.bond_weights()) {
      typename WEIGHT_TABLE::bond_weight_iterator bwi, bwi_end;
      for (boost::tie(bwi, bwi_end) = wt.bond_weights(); bwi != bwi_end; ++bwi) {
        if (bwi->first >= sharing.bond_offset(tid) &&
            bwi->first < sharing.bond_offset(tid) + sharing.num_bonds_local(tid)) {
          for (int g = 0; g < WEIGHT_TABLE::num_bond_graphs; ++g) {
            if (is_nonzero<1>(bwi->second.v[g])) {
              if (bond_graph_t::is_valid_gid(g)) {
                graph_[tid].push_back(local_graph_t(g, location_t::bond_location(bwi->first)));
                w[tid].push_back(bwi->second.v[g]);
                weight_[tid] += bwi->second.v[g];
              } else {
                boost::throw_exception(std::runtime_error("graph_chooser::init"));
              }
            }
          }
        }
      }
      if (w[tid].size()) dist_graph_[tid].init(w[tid]);
      if (is_zero<2>(weight_[tid])) weight_[tid] = 1.0e-20;
    } // end omp parallel
    //// BOOST_FOREACH(typename WEIGHT_TABLE::bond_weight_t const& bw, wt.bond_weights()) {
    typename WEIGHT_TABLE::bond_weight_iterator bwi, bwi_end;
    for (boost::tie(bwi, bwi_end) = wt.bond_weights(); bwi != bwi_end; ++bwi) {
      if (!is_path_integral && bond_graph_t::has_nontrivial_diagonal_choice_helper())
        diag_.push_back(typename bond_graph_t::diagonal_choice_helper(bwi->second.v));
      if (bond_graph_t::has_nontrivial_offdiagonal_choice_helper())
        offdiag_.push_back(typename bond_graph_t::offdiagonal_choice_helper(bwi->second.v));
    }
  }
#endif

  template<typename RNG>
  const local_graph_t& graph(RNG& rng) const { return graph_[0][dist_graph_[0](rng)]; }
  template<typename RNG>
  const local_graph_t& graph(RNG& rng, int tid) const { return graph_[tid][dist_graph_[tid](rng)]; }
  template<typename RNG>
  local_graph_t diagonal(RNG& /* rng */, const location_t& loc, int /* c */) const {
    return local_graph_t(0, loc);
  }
  template<typename RNG>
  local_graph_t diagonal(RNG& rng, const location_t& loc, int c0, int c1) const {
    return local_graph_t(bond_graph_t::choose_diagonal(rng, diag_, pos(loc), c0, c1), loc);
  }
  template<typename RNG>
  local_graph_t offdiagonal(RNG& rng, const location_t& loc, int c0, int c1) const {
    return local_graph_t(bond_graph_t::choose_offdiagonal(rng, offdiag_, pos(loc), c0, c1), loc);
  }
  double weight() const { return weight_[0]; }
  double weight(int tid) const { return weight_[tid]; }

private:
  std::vector<random_choice_walker_d<> > dist_graph_;
  std::vector<double> weight_;
  std::vector<std::vector<local_graph_t> > graph_;
  std::vector<typename bond_graph_t::diagonal_choice_helper> diag_;
  std::vector<typename bond_graph_t::offdiagonal_choice_helper> offdiag_;
};

} // end namespace looper

#endif // LOOPER_GARPH_IMPL_H
