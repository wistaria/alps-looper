/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_PATH_INTEGRAL_H
#define LOOPER_PATH_INTEGRAL_H

#include <looper/amida.h>
#include <looper/fill_duration.h>
#include <looper/node.h>
#include <looper/permutation.h>
#include <looper/union_find.h>
#include <looper/virtual_graph.h>
#include <looper/weight.h>

#include <boost/integer_traits.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/throw_exception.hpp>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace looper {

template<bool HasCTime = false> class pi_node;

template<>
class pi_node<false> : public qmc_node
{
public:
  typedef qmc_node                base_type;
  typedef double                  time_type;
  typedef std::complex<time_type> ctime_type;
  BOOST_STATIC_CONSTANT(bool, has_ctime = false);

  pi_node() : base_type(), time_(0.) {}

  time_type time() const { return time_; }
  ctime_type ctime() const { return std::exp(2 * M_PI * time_); }
  void set_time(time_type t) { time_ = t; }

  std::ostream& output(std::ostream& os) const {
    return base_type::output(os) << " time = " << time_;
  }
  alps::ODump& save(alps::ODump& od) const {
    return base_type::save(od) << time_;
  }
  alps::IDump& load(alps::IDump& id) {
    return base_type::load(id) >> time_;
  }

private:
  time_type time_;
};

template<>
class pi_node<true> : public pi_node<false>
{
public:
  typedef pi_node<false>          base_type;
  typedef base_type::time_type    time_type;
  typedef std::complex<time_type> ctime_type;
  BOOST_STATIC_CONSTANT(bool, has_ctime = true);

  pi_node() : base_type(), ctime_() {}

  ctime_type ctime() const { return ctime_; }
  void set_time(time_type t) {
    base_type::set_time(t);
    ctime_ = std::exp(2 * M_PI * t);
  }

  std::ostream& output(std::ostream& os) const {
    return base_type::output(os) << " ctime = " << ctime_;
  }
  alps::ODump& save(alps::ODump& od) const { return base_type::save(od); }
  alps::IDump& load(alps::IDump& id) {
    base_type::load(id);
    ctime_ = std::exp(2 * M_PI * time());
    return id;
  }

private:
  ctime_type ctime_;
};


template<class G, class M, class W = weight::xxz, class N = pi_node<> >
struct path_integral;

template<class G, class M, class W, class N>
struct path_integral<virtual_graph<G>, M, W, N>
{
  typedef virtual_graph<G>                      vg_type;
  typedef typename virtual_graph<G>::graph_type graph_type;
  typedef M                                     model_type;
  typedef W                                     weight_type;
  typedef N                                     node_type;

  typedef typename boost::graph_traits<graph_type>::edge_iterator
                                                    edge_iterator;
  typedef typename boost::graph_traits<graph_type>::edge_descriptor
                                                    edge_descriptor;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
                                                    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::vertex_descriptor
                                                    vertex_descriptor;

  struct parameter_type
  {
    typedef path_integral<virtual_graph<G>, M, W, N> qmc_type;

    typedef virtual_graph<G>               vg_type;
    typedef typename vg_type::graph_type   graph_type;
    typedef typename vg_type::mapping_type mapping_type;
    typedef M                              model_type;
    typedef W                              weight_type;

    template<class RG>
    parameter_type(const RG& rg, const model_type& m, double b, double fs)
      : virtual_graph(), model(m), beta(b), is_bipartite(false), weight()
    {
      generate_virtual_graph(virtual_graph, rg, model);
      is_bipartite = alps::set_parity(virtual_graph.graph);
      
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(virtual_graph.graph);
           ei != ei_end; ++ei) {
        weight.push_back(weight_type(
          model.bond(bond_type(*ei, virtual_graph.graph)), fs));
      }
    }

    vg_type                  virtual_graph;
    const model_type&        model;
    double                   beta;
    bool                     is_bipartite;
    std::vector<weight_type> weight;
  };

  struct config_type
  {
    typedef path_integral<virtual_graph<G>, M, W, N> qmc_type;

    typedef N                                node_type;
    typedef amida<node_type>                 wl_type;
    typedef typename wl_type::iterator       iterator;
    typedef typename wl_type::const_iterator const_iterator;

    BOOST_STATIC_CONSTANT(bool, is_path_integral = true);
    BOOST_STATIC_CONSTANT(bool, is_sse = false);

    wl_type wl;
    unsigned int num_loops0;
    unsigned int num_loops;

    alps::ODump& save(alps::ODump& od) const { return od << wl; }
    alps::IDump& load(alps::IDump& id) { return id >> wl; }
  };

  // helper functions: segment_d, segment_u

  static typename node_type::segment_type&
  segment_d(const typename config_type::iterator& itr)
  {
    if (itr.at_boundary()) {
      return itr->loop_segment(0);
    } else {
      if (itr->is_refl()) {
        return itr->loop_segment(1);
      } else {
        return itr->loop_segment(itr.leg());
      }
    }
  }

  static const typename node_type::segment_type&
  segment_d(const typename config_type::const_iterator& itr)
  {
    if (itr.at_boundary()) {
      return itr->loop_segment(0);
    } else {
      if (itr->is_refl()) {
        return itr->loop_segment(1);
      } else {
        return itr->loop_segment(itr.leg());
      }
    }
  }

  static typename node_type::segment_type&
  segment_u(const typename config_type::iterator& itr)
  {
    if (itr.at_boundary()) {
      return itr->loop_segment(0);
    } else {
      if (itr->is_refl()) {
        return itr->loop_segment(0);
      } else {
        return itr->loop_segment(1-itr.leg());
      }
    }
  }

  static const typename node_type::segment_type&
  segment_u(const typename config_type::const_iterator& itr)
  {
    if (itr.at_boundary()) {
      return itr->loop_segment(0);
    } else {
      if (itr->is_refl()) {
        return itr->loop_segment(0);
      } else {
        return itr->loop_segment(1-itr.leg());
      }
    }
  }

  //
  // update functions
  //

  // initialize
  static void initialize(config_type& config, const vg_type& vg)
  {
    config.wl.init(boost::num_vertices(vg.graph));

    vertex_iterator vi_end = boost::vertices(vg.graph).second;
    for (vertex_iterator vi = boost::vertices(vg.graph).first;
         vi != vi_end; ++vi) {
      // all up
      config.wl.series(*vi).first ->conf() = 0;
      config.wl.series(*vi).second->conf() = 0;
      config.wl.series(*vi).first ->set_time(0.);
      config.wl.series(*vi).second->set_time(1.);
    }
  }

  static void initialize(config_type& config,
                         const parameter_type& p)
  { initialize(config, p.virtual_graph); }

  template<class RNG>
  static void generate_loops(config_type& config, const vg_type& vg,
                             const model_type& model, double beta,
                             const std::vector<weight_type>& weight,
                             RNG& uniform_01)
  {
    typedef typename config_type::iterator  iterator;
    typedef typename config_type::node_type node_type;

    //
    // labeling
    //

    {
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(vg.graph);
           ei != ei_end; ++ei) {
        int bond = boost::get(edge_index_t(), vg.graph, *ei); // bond index

        // setup iterators
        iterator itr0 = config.wl.series(boost::source(*ei, vg.graph)).first;
        iterator itr1 = config.wl.series(boost::target(*ei, vg.graph)).first;
        int c0 = itr0->conf();
        int c1 = itr1->conf();
        ++itr0;
        ++itr1;

        std::vector<double> trials;
        fill_duration(uniform_01, trials, beta * weight[bond].density());

        // iteration up to t = 1
        std::vector<double>::const_iterator ti_end = trials.end();
        for (std::vector<double>::const_iterator ti = trials.begin();
             ti != ti_end; ++ti) {
          while (itr0->time() < *ti) {
            if (itr0->bond() == bond) // labeling existing link
              itr0->set_old(uniform_01() < weight[bond].p_reflect());
            if (itr0->is_old()) c0 ^= 1;
            ++itr0;
          }
          while (itr1->time() < *ti) {
            if (itr1->is_old()) c1 ^= 1;
            ++itr1;
          }
          if (uniform_01() < weight[bond].p_accept(c0, c1)) {
            // insert new link
            iterator itr_new =
              config.wl.insert_link_prev(node_type(), itr0, itr1).first;
            itr_new->set_time(*ti);
            itr_new->set_bond(bond);
            itr_new->set_new((c0 ^ c1), 
                             (uniform_01() < weight[bond].p_freeze(c0 ^ c1)));
          }
        }
        while (!itr0.at_top()) {
          if (itr0->bond() == bond)        // labeling existing link
            itr0->set_old(uniform_01() < weight[bond].p_reflect());
          ++itr0;
        }
      }
    }

    //
    // cluster identification using union-find algorithm
    //

    {
      vertex_iterator vi_end = boost::vertices(vg.graph).second;
      for (vertex_iterator vi = boost::vertices(vg.graph).first;
           vi != vi_end; ++vi) {
        // setup iterators
        iterator itrD = config.wl.series(*vi).first;
        iterator itrU = boost::next(itrD);

        // iteration up to t = 1
        for (;; itrD = itrU, ++itrU) {
          // Intel C++ (icc) 8.0 does not understand (;; itrD = itrU++)
          // correctly when -O1 or higher is specified.  -- ST 2004.04.01
          union_find::unify(segment_d(itrU), segment_u(itrD));
          if (itrU.at_top()) break; // finish
          if (itrU.leg() == 0 && itrU->is_frozen()) // frozen link
            union_find::unify(itrU->loop_segment(0), itrU->loop_segment(1));
        }
      }
    }

    // connect bottom and top with random permutation
    {
      std::vector<int> r, c0, c1;
      for (int i = 0; i < vg.mapping.num_groups(); ++i) {
        int s2 = vg.mapping.num_virtual_vertices(i);
        int offset = *(vg.mapping.virtual_vertices(i).first);
        if (s2 == 1) {
          // S=1/2: just connect top and bottom
          union_find::unify(config.wl.series(offset).first ->loop_segment(0),
                            config.wl.series(offset).second->loop_segment(0));
        } else {
          // S>=1
          r.resize(s2);
          c0.resize(s2);
          c1.resize(s2);
          vertex_iterator vi_end = vg.mapping.virtual_vertices(i).second;
          for (vertex_iterator vi = vg.mapping.virtual_vertices(i).first;
               vi != vi_end; ++vi) {
            r[*vi - offset] = *vi - offset;
            c0[*vi - offset] = config.wl.series(*vi).first->conf();
            c1[*vi - offset] = config.wl.series(*vi).second->conf();
          }
          restricted_random_shuffle(r.begin(), r.end(),
                                    c0.begin(), c0.end(),
                                    c1.begin(), c1.end(),
                                    uniform_01);
          for (int j = 0; j < s2; ++j) {
            union_find::unify(
              config.wl.series(offset +   j ).first ->loop_segment(0),
              config.wl.series(offset + r[j]).second->loop_segment(0));
          }
        }
      }
    }

    // counting and indexing loops
    {
      config.num_loops0 = 0;
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph);
           vi != vi_end; ++vi) {
        iterator itrB, itrT;
        boost::tie(itrB, itrT) = config.wl.series(*vi);
        if (itrB->loop_segment(0).index == loop_segment::undefined) {
          if (itrB->loop_segment(0).root()->index == loop_segment::undefined)
            itrB->loop_segment(0).root()->index = (config.num_loops0)++;
          itrB->loop_segment(0).index = itrB->loop_segment(0).root()->index;
        }
        if (itrT->loop_segment(0).index == loop_segment::undefined) {
          if (itrT->loop_segment(0).root()->index == loop_segment::undefined)
            itrT->loop_segment(0).root()->index = (config.num_loops0)++;
          itrT->loop_segment(0).index = itrT->loop_segment(0).root()->index;
        }
      }
    }

    {
      config.num_loops = config.num_loops0;
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph);
           vi != vi_end; ++vi) {
        // setup iterator
        iterator itr = boost::next(config.wl.series(*vi).first);

        // iteration up to t = 1
        while (!itr.at_top()) {
          if (itr.leg() == 0) {
            if (itr->loop_segment(0).index == loop_segment::undefined) {
              if (itr->loop_segment(0).root()->index ==
                  loop_segment::undefined)
                itr->loop_segment(0).root()->index = (config.num_loops)++;
              itr->loop_segment(0).index = itr->loop_segment(0).root()->index;
            }
            if (itr->loop_segment(1).index == loop_segment::undefined) {
              if (itr->loop_segment(1).root()->index ==
                  loop_segment::undefined)
                itr->loop_segment(1).root()->index = (config.num_loops)++;
              itr->loop_segment(1).index = itr->loop_segment(1).root()->index;
            }
          }
          ++itr;
        }
      }
    }
  }

  template<class RNG>
  static void generate_loops(config_type& config,
                             const parameter_type& p,
                             RNG& uniform_01)
  {
    generate_loops(config, p.virtual_graph, p.model, p.beta, p.weight,
                   uniform_01);
  }

  template<class RNG>
  static void flip_and_cleanup(config_type& config,
                               const vg_type& vg, RNG& uniform_01)
  {
    typedef typename config_type::iterator iterator;
    typedef RNG rng_type;

    std::vector<int> flip(config.num_loops);
    {
      boost::variate_generator<rng_type&, boost::uniform_smallint<> >
        uniform_int01(uniform_01, boost::uniform_smallint<>(0, 1));
      std::vector<int>::iterator itr = flip.begin();
      std::vector<int>::iterator itr_end = flip.end();
      for (; itr != itr_end; ++itr) *itr = uniform_int01();
    }
    // Intel C++ (icc) 8.0 does not understand the following lines
    // correctly when -xW vectorized option is specified.  -- ST 2004.04.01
    //
    // std::generate(flip.begin(), flip.end(),
    //   boost::variate_generator<rng_type&, boost::uniform_smallint<> >(
    //   uniform_01, boost::uniform_smallint<>(0, 1)));

    // flip spins
    {
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph);
           vi != vi_end; ++vi) {
        iterator itrB, itrT;
        boost::tie(itrB, itrT) = config.wl.series(*vi);
        if (flip[itrB->loop_segment(0).index] == 1) itrB->flip_conf();
        itrB->clear_graph();
        if (flip[itrT->loop_segment(0).index] == 1) itrT->flip_conf();
        itrT->clear_graph();
      }
    }

    // upating links
    {
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(vg.graph);
           vi != vi_end; ++vi) {
        // setup iterator
        iterator itr = boost::next(config.wl.series(*vi).first);

        // iteration up to t = 1
        while (!itr.at_top()) {
          if (itr.leg() == 0) {
            if (itr->is_new() ^
                (flip[itr->loop_segment(0).index] ^
                 flip[itr->loop_segment(1).index] == 1)) {
              ++itr;
              config.wl.erase(boost::prior(itr));
            } else {
              itr->clear_graph();
              ++itr;
            }
          } else {
            ++itr;
          }
        }
      }
    }
  }
  template<class RNG>
  static void flip_and_cleanup(config_type& config,
                               const parameter_type& p,
                               RNG& uniform_01)
  { flip_and_cleanup(config, p.virtual_graph, uniform_01); }

  // measurements

  static int loop_index_0(int i, const config_type& config)
  { return config.wl.series(i).first->loop_segment(0).index; }

  static int loop_index_1(int i, const config_type& config)
  { return config.wl.series(i).second->loop_segment(0).index; }

  static double static_sz(int i, const config_type& config)
  { return 0.5 - double(config.wl.series(i).first->conf()); }

  static double dynamic_sz(int i, const config_type& config, const vg_type&)
  {
    typedef typename config_type::const_iterator const_iterator;

    double sz = 0.;

    // setup iterators
    const_iterator itrD = config.wl.series(i).first;
    const_iterator itrU = boost::next(itrD);
    double c = static_sz(i, config);

    // iteration up to t = 1
    for (;; itrD = itrU, ++itrU, c *= -1) {
      sz += c * (itrU->time() - itrD->time());
      if (itrU.at_top()) break;
    }
    return sz;
  }

  static double static_correlation(const config_type& config,
                                   const vertex_descriptor& v0,
                                   const vertex_descriptor& v1)
  {
    typedef typename config_type::const_iterator const_iterator;

    double corr = 0.;
    double t = 0.;

    // setup iterators
    const_iterator itr0 = config.wl.series(v0).first;
    const_iterator itr1 = config.wl.series(v1).first;
    double c = static_sz(v0, config) * static_sz(v1, config);
    while (true) {
      if (itr0->time() < itr1->time()) {
        corr += c * (itr0->time() - t);
        if (itr0.at_top()) break;
        t = itr0->time();
        ++itr0;
      } else {
        corr += c * (itr1->time() - t);
        if (itr1.at_top()) break;
        t = itr1->time();
        ++itr1;
      }
      c *= -1.;
    }
    return corr;
  }

  static int num_offdiagonals(const config_type& config)
  { return config.wl.num_links(); }

  static double energy_offset(const vg_type& vg, const M& model)
  {
    double offset = 0.;
    edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(vg.graph); ei != ei_end; ++ei)
      offset += model.bond(bond_type(*ei, vg.graph)).c();
    return offset / vg.num_real_edges;
  }
  static double energy_offset(const parameter_type& p)
  { return energy_offset(p.virtual_graph, p.model); }

}; // struct path_integral

} // namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<bool HasCTime>
std::ostream& operator<<(std::ostream& os,
                         const looper::pi_node<HasCTime>& n) {
  n.output(os);
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace looper
#endif

#endif // LOOPER_PATH_INTEGRAL_H
