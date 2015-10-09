/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
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

// Weighted Union-Find Algorithm
// Reference:
//   D. Knuth,
//   `The Art of Computer Programming, Vol. 1, Fundamental Algorithms'
//   3rd edition (Addison Wesley, Reading, 1997) Sec 2.3.3.

#ifndef LOOPER_UNION_FIND_H
#define LOOPER_UNION_FIND_H

#ifndef ALPS_INDEP_SOURCE
# include <alps/config.h>
# if defined(LOOPER_ENABLE_OPENMP) && defined(ALPS_ENABLE_OPENMP_WORKER) && !defined(LOOPER_OPENMP)
#  define LOOPER_OPENMP
# endif
#else
# if defined(LOOPER_ENABLE_OPENMP) && defined(_OPENMP) && !defined(LOOPER_OPENMP)
#  define LOOPER_OPENMP
# endif
#endif

#include <algorithm> // for std::swap
#include <vector>
#include <iostream>

#ifdef LOOPER_OPENMP
# include <omp.h>
# include "atomic.h"
#endif

namespace looper {
namespace union_find {

class node {
public:
  node() : parent_(-1) {} // root node with weight = 1
  ~node() {}
  bool is_root() const { return parent_ <= 0; }
  void set_parent(int parent) { parent_ = parent + 1; }
  int parent() const { return parent_ - 1; }
  void set_weight(int w) { parent_ = -w; }
  int weight() const { return -parent_; }
  void set_id(int id) { id_ = id; }
  int id() const { return id_; }
#ifdef LOOPER_OPENMP
  int lock_root() {
    int p = parent_;
    if (p < 0 && compare_and_swap(parent_, p, 0)) {
      return -p;
    } else {
      return 0;
    }
  }
  // unlock can be done by set_parent or set_weight
#endif
private:
  int parent_; // negative for root fragment
  int id_;
};

class node_c {
public:
  node_c() : parent_(1 ^ id_mask) {} // root node with weight = 1
  ~node_c() {}
  int is_root() const { return (int) parent_ <= 0; }
  void set_parent(int parent) { parent_ = parent + 1; }
  int parent() const { return parent_ - 1; }
  void set_weight(int w) { parent_ = w ^ id_mask; }
  int weight() const { return parent_ ^ id_mask; }
  void set_id(int id) { parent_ = id ^ id_mask; }
  int id() const { return parent_ ^ id_mask; }
#ifdef LOOPER_OPENMP
  int lock_root() {
    int p = parent_;
    if (p < 0 && compare_and_swap(parent_, p, 0)) {
      return -p;
    } else {
      return 0;
    }
  }
  // unlock can be done by set_parent or set_weight
#endif
private:
  static const int id_mask = -1;
  int parent_; // non-positive for root fragment
};

class node_noweight {
public:
  node_noweight() : parent_(id_mask) {} // root note with id = 0
  ~node_noweight() {}
  int is_root() const { return (int) parent_ <= 0; }
  void set_parent(int parent) { parent_ = parent + 1; }
  int parent() const { return parent_ - 1; }
  void set_weight(int) { set_id(0); } // dummy routine for unlock
  int weight() const { return 0; } // dummy
  void set_id(int id) { parent_ = id ^ id_mask; }
  int id() const { return parent_ ^ id_mask; }
#ifdef LOOPER_OPENMP
  int lock_root() {
    int p = parent_;
    if (p < 0 && compare_and_swap(parent_, p, 0)) {
      return 1;
    } else {
      return 0;
    }
  }
  // unlock can be done by set_parent, set_weight or set_id
#endif
private:
  static const int id_mask = -1;
  int parent_; // negative for root fragment
};

// thread-unsafe
template<class T>
inline int add(std::vector<T>& v) {
  v.push_back(T());
  return v.size() - 1; // return index of new node
}

template<class T>
inline int root_index(std::vector<T> const& v, int g) {
#ifdef LOOPER_OPENMP
  T c = v[g];
  while (!c.is_root()) {
    g = c.parent();
    c = v[g];
  }
#else
  while (!v[g].is_root()) {
    g = v[g].parent();
  }
#endif
  return g;
}

// root_index with path-halving
// Note: this is not thread-safe, but is really safe as long as called from unify_*
template<class T>
inline int root_index_ph(std::vector<T>& v, int g) {
  if (v[g].is_root()) return g;
  while (true) {
    int p = v[g].parent();
    if (v[p].is_root()) return p;
    v[g].set_parent(v[p].parent());
    g = p;
  }
}

template<class T>
inline T const& root(std::vector<T> const& v, int g) { return v[root_index(v, g)]; }

template<class T>
inline T const& root(std::vector<T> const& v, T const& n) {
  return n.is_root() ? n : root(v, n.parent());
}

template<class T>
inline int cluster_id(std::vector<T> const& v, int g) { return root(v, g).id(); }

template<class T>
inline int cluster_id(std::vector<T> const& v, T const& n) { return root(v, n).id(); }

template<class T>
void set_root(std::vector<T>& v, int g) {
#ifdef LOOPER_OPENMP
  while(true) {
    int r = root_index(v, g);
    if (r == g) {
      return;
    } else {
      int w = v[r].lock_root();
      if (w != 0) {
        v[g].set_weight(w);
        v[r].set_parent(g); // release lock
        return;
      }
    }
  }
#else
  int r = root_index(v, g);
  if (r != g) {
    v[g] = v[r];
    v[r].set_parent(g);
  }
#endif
}

template<class T>
inline void update_link(std::vector<T>& v, int g, int r) {
  while (g != r) {
    int p = v[g].parent();
    v[g].set_parent(r);
    g = p;
  }
}

// WARNING: this is not thread-safe
template<class T>
inline int unify_compress(std::vector<T>& v, int g0, int g1) {
  using std::swap;
  int r0 = root_index(v, g0);
  int r1 = root_index(v, g1);
  if (r0 != r1) {
#if !defined(LOOPER_USE_DETERMINISTIC_UNIFY)
    if (v[r0].weight() < v[r1].weight()) swap(r0, r1);
#else
    if (r0 > r1) swap(r0, r1);
#endif
    v[r0].set_weight(v[r0].weight() + v[r1].weight());
    v[r1].set_parent(r0);
  }
  update_link(v, g0, r0);
  update_link(v, g1, r0);
  return r0; // return (new) root node
}

template<class T>
inline int unify_pathhalving(std::vector<T>& v, int g0, int g1) {
  using std::swap;
  int r0 = root_index_ph(v, g0);
  int r1 = root_index_ph(v, g1);
#ifdef LOOPER_OPENMP
  int w0 = 0;
  int w1 = 0;
  while (true) {
    if (r1 == r0) return r0; // g0 and g1 belong to the same cluster
    w0 = v[r0].lock_root();
    w1 = v[r1].lock_root();
    if (w0 != 0 && w1 != 0) break;
    if (w0 != 0) v[r0].set_weight(w0); // release lock
    if (w1 != 0) v[r1].set_weight(w1); // release lock
    r0 = root_index_ph(v, r0);
    r1 = root_index_ph(v, r1);
  }
# if !defined(LOOPER_USE_DETERMINISTIC_UNIFY)
  if (w0 < w1) swap(r0, r1);
# else
  if (r0 > r1) swap(r0, r1);
# endif
  v[r0].set_weight(w0 + w1); // release lock
  v[r1].set_parent(r0); // release lock
#else
  if (r0 != r1) {
# if !defined(LOOPER_USE_DETERMINISTIC_UNIFY)
    if (v[r0].weight() < v[r1].weight()) swap(r0, r1);
# else
    if (r0 > r1) swap(r0, r1);
# endif
    v[r0].set_weight(v[r0].weight() + v[r1].weight());
    v[r1].set_parent(r0);
  }
#endif
  return r0; // return (new) root node
}

template<class T>
inline int unify(std::vector<T>& v, int g0, int g1) {
  return unify_pathhalving(v, g0, g1);
}

template<class T>
inline void output(std::vector<T> const& v, std::ostream& os = std::cout) {
  for (int i = 0; i < v.size(); ++i) {
    os << "node " << i << ": ";
    int g = i;
    if (v[g].is_root()) {
      os << "root (id = " << v[g].id() << ")" << std::endl;
    } else {
      while (true) {
        if (v[g].is_root()) {
          os << g << std::endl;
          break;
        } else {
          os << g << " -> ";
          g = v[g].parent();
        }
      }
    }
  }
}

template<typename T>
int count_root(std::vector<T>& v, int start, int n) {
  int nc = 0;
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) ++nc;
  return nc;
}

template<typename T>
int count_root_p(std::vector<T>& v, int start, int n) {
  int nc = 0;
  #pragma omp for schedule(static) nowait
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) ++nc;
  return nc;
}

template<typename T>
int set_id(std::vector<T>& v, int start, int n, int nc) {
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) v[i].set_id(nc++);
  return nc;
}

template<typename T>
int set_id_p(std::vector<T>& v, int start, int n, int nc) {
  #pragma omp for schedule(static) nowait
  for (int i = start; i < start + n; ++i)
    if (v[i].is_root()) v[i].set_id(nc++);
  return nc;
}

template<typename T>
void copy_id(std::vector<T>& v, int start, int n) {
  for (int i = start; i < start + n; ++i)
    v[i].set_id(cluster_id(v, i));
}

template<typename T>
void copy_id_p(std::vector<T>& v, int start, int n) {
  #pragma omp for schedule(static) nowait
  for (int i = start; i < start + n; ++i)
    v[i].set_id(cluster_id(v, i));
}

template<typename T>
inline void pack_tree(std::vector<T>& v, int n) {
#ifdef LOOPER_OPENMP
  int g, w; // workaround for FCC OpenMP bug? -- ST 2010-11-22
  #pragma omp parallel for schedule(static) private(g, w)
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      g = v[i].parent();
      while (true) {
        if (g < n) {
          // encounter node with index < n
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n
          w = v[g].lock_root();
          if (w != 0) {
            v[i].set_weight(w);
            v[g].set_parent(i); // release lock
            break;
          }
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      int g = v[i].parent();
      while (true) {
        if (g < n) {
          // encounter node with index < n
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n
          v[i] = v[g];
          v[g].set_parent(i);
          break;
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#endif
}

// pack tree so that nodes with id [0...n) and [m...) come upper
template<typename T>
inline void pack_tree(std::vector<T>& v, int n, int m) {
#ifdef LOOPER_OPENMP
  int g, w; // workaround for FCC OpenMP bug? -- ST 2010-11-22
  #pragma omp parallel for schedule(static) private(g, w)
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          w = v[g].lock_root();
          if (w != 0) {
            v[i].set_weight(w);
            v[g].set_parent(i); // release lock
            break;
          }
        } else {
          g = v[g].parent();
        }
      }
    }
  }
  #pragma omp parallel for schedule(static) private(g, w)
  for (int i = m; i < m + n; ++i) {
    if (!v[i].is_root()) {
      g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          w = v[g].lock_root();
          if (w != 0) {
            v[i].set_weight(w);
            v[g].set_parent(i); // release lock
            break;
          }
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    if (!v[i].is_root()) {
      int g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          v[i] = v[g];
          v[g].set_parent(i);
          break;
        } else {
          g = v[g].parent();
        }
      }
    }
  }
  for (int i = m; i < m + n; ++i) {
    if (!v[i].is_root()) {
      int g = v[i].parent();
      while (true) {
        if (g < n || g >= m) {
          // encounter node with index < n or >= m
          v[i].set_parent(g);
          break;
        } else if (v[g].is_root()) {
          // found root with index >= n and < m
          v[i] = v[g];
          v[g].set_parent(i);
          break;
        } else {
          g = v[g].parent();
        }
      }
    }
  }
#endif
}

} // end namespace union_find
} // end namespace looper

#endif
