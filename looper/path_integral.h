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
#include <looper/graph.h>
#include <looper/node.h>
#include <looper/permutation.h>
#include <looper/sign.h>
#include <looper/union_find.h>
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
struct path_integral
{
  typedef G graph_type;
  typedef M model_type;
  typedef W weight_type;
  typedef N node_type;

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
    typedef path_integral<G, M, W, N> qmc_type;

    typedef G                           graph_type;
    typedef virtual_mapping<graph_type> mapping_type;
    typedef M                           model_type;
    typedef W                           weight_type;

    template<class RG>
    parameter_type(const alps::graph_helper<RG>& gh, const model_type& m,
                   double b, double fs)
      : model(m), vgraph(), vmap(), 
        num_real_vertices(boost::num_vertices(rg)),
        num_real_edges(boost::num_edges(rg)),
        beta(b), sz_conserved(true), is_bipartite(false), weight(),
        ez_offset(0.0)
    {
      generate_virtual_graph(rg, model, vgraph, vmap);
      is_bipartite = alps::set_parity(vgraph);

      // if (model.is_signed() || model.is_classically_frustrated())
      //   fs = std::max(fs, 0.1);
      if (model.is_classically_frustrated()) fs = std::max(fs, 0.1);

      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(vgraph); ei != ei_end; ++ei) {
        weight.push_back(weight_type(model.bond(bond_type(*ei, vgraph)), fs));
        ez_offset += model.bond(bond_type(*ei, vgraph)).c();
      }
    }

    const model_type&        model;
    graph_type               vgraph;
    mapping_type             vmap;
    unsigned int             num_real_vertices;
    unsigned int             num_real_edges;
    double                   beta;
    bool                     sz_conserved;
    bool                     is_bipartite;
    std::vector<weight_type> weight;
    double                   ez_offset;
  };

  struct config_type : public sign_info
  {
    typedef path_integral<G, M, W, N> qmc_type;

    typedef N                                node_type;
    typedef amida<node_type>                 wl_type;
    typedef typename wl_type::iterator       iterator;
    typedef typename wl_type::const_iterator const_iterator;

    BOOST_STATIC_CONSTANT(bool, is_path_integral = true);
    BOOST_STATIC_CONSTANT(bool, is_sse = false);

    config_type() : sign_info() {}

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
  static void initialize(config_type& config, const graph_type& vgraph)
  {
    config.wl.init(boost::num_vertices(vgraph));

    vertex_iterator vi_end = boost::vertices(vgraph).second;
    for (vertex_iterator vi = boost::vertices(vgraph).first;
         vi != vi_end; ++vi) {
      // all up
      config.wl.series(*vi).first ->conf() = 0;
      config.wl.series(*vi).second->conf() = 0;
      config.wl.series(*vi).first ->set_time(0.);
      config.wl.series(*vi).second->set_time(1.);
    }

    config.sign = 1;
  }

  static void initialize(config_type& config, const parameter_type& param)
  { initialize(config, param.vgraph); }

  template<class RNG>
  static void generate_loops(config_type& config,
                             const parameter_type& param,
                             RNG& uniform_01)
  {
    typedef typename config_type::iterator  iterator;
    typedef typename config_type::node_type node_type;

    //
    // labeling
    //

    {
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(param.vgraph);
           ei != ei_end; ++ei) {
        int bond = boost::get(boost::edge_index, param.vgraph, *ei);

        // setup iterators
        iterator itr0 =
          config.wl.series(boost::source(*ei, param.vgraph)).first;
        iterator itr1 =
          config.wl.series(boost::target(*ei, param.vgraph)).first;
        int c0 = itr0->conf();
        int c1 = itr1->conf();
        ++itr0;
        ++itr1;

        std::vector<double> trials;
        fill_duration(uniform_01, trials, param.beta * param.weight[bond].density());

        // iteration up to t = 1
        std::vector<double>::const_iterator ti_end = trials.end();
        for (std::vector<double>::const_iterator ti = trials.begin();
             ti != ti_end; ++ti) {
          while (itr0->time() < *ti) {
            if (itr0->bond() == bond) // labeling existing link
              itr0->set_old(uniform_01() < param.weight[bond].p_reflect());
            if (itr0->is_old()) c0 ^= 1;
            ++itr0;
          }
          while (itr1->time() < *ti) {
            if (itr1->is_old()) c1 ^= 1;
            ++itr1;
          }
          if (uniform_01() < param.weight[bond].p_accept(c0, c1)) {
            // insert new link
            iterator itr_new =
              config.wl.insert_link_prev(node_type(), itr0, itr1).first;
            itr_new->set_time(*ti);
            itr_new->set_bond(bond);
            itr_new->set_new((c0 ^ c1), 
                             (uniform_01() < param.weight[bond].p_freeze(c0 ^ c1)));
          }
        }
        while (!itr0.at_top()) {
          if (itr0->bond() == bond)        // labeling existing link
            itr0->set_old(uniform_01() < param.weight[bond].p_reflect());
          ++itr0;
        }
      }
    }

    //
    // cluster identification using union-find algorithm
    //

    {
      vertex_iterator vi_end = boost::vertices(param.vgraph).second;
      for (vertex_iterator vi = boost::vertices(param.vgraph).first;
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
      for (int i = 0; i < param.vmap.num_groups(); ++i) {
        int s2 = param.vmap.num_virtual_vertices(i);
        int offset = *(param.vmap.virtual_vertices(i).first);
        if (s2 == 1) {
          // S=1/2: just connect top and bottom
          union_find::unify(config.wl.series(offset).first ->loop_segment(0),
                            config.wl.series(offset).second->loop_segment(0));
        } else {
          // S>=1
          r.resize(s2);
          c0.resize(s2);
          c1.resize(s2);
          vertex_iterator vi_end = param.vmap.virtual_vertices(i).second;
          for (vertex_iterator vi = param.vmap.virtual_vertices(i).first;
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
      for (boost::tie(vi, vi_end) = boost::vertices(param.vgraph);
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
      for (boost::tie(vi, vi_end) = boost::vertices(param.vgraph);
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

    // negative sign information
    if (param.model.is_signed()) {
      config.num_merons0 = 0;
      config.num_merons = 0;
      config.loop_sign.clear();
      config.loop_sign.resize(config.num_loops, 0);
      edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(param.vgraph);
           ei != ei_end; ++ei) {
        int bond = boost::get(boost::edge_index, param.vgraph,
                              *ei);
        if (param.weight[bond].sign() < 0) {
          typename config_type::const_iterator
            itr = boost::next(config.wl.series(
              boost::source(*ei, param.vgraph)).first);
          while (!itr.at_top()) {
            if (itr->bond() == bond) {
              ++config.loop_sign[itr->loop_segment(0).index];
              ++config.loop_sign[itr->loop_segment(1).index];
            }
            ++itr;
          }
        }
      }
      for (int i = 0; i < config.num_loops0; ++i) {
        // sign = 1 for even number of negative links
        //       -1 for odd
        if (config.loop_sign[i] % 2 == 1) {
          config.loop_sign[i] = -1;
          ++config.num_merons0;
          ++config.num_merons;
        } else {
          config.loop_sign[i] = 1;
        }
      }
      for (int i = config.num_loops0; i < config.num_loops; ++i) {
        // sign = 1 for even number of negative links
        //       -1 for odd
        if (config.loop_sign[i] % 2 == 1) {
          config.loop_sign[i] = -1;
          ++config.num_merons;
        } else {
          config.loop_sign[i] = 1;
        }
      }
    }
  }

  template<class RNG>
  static void flip_and_cleanup(config_type& config,
                               const parameter_type& param,
                               RNG& uniform_01)
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
    // Intel C++ (icc) 8.x does not understand the following lines
    // correctly when -xW vectorized option is specified.  -- ST 2004.04.01
    //
    // std::generate(flip.begin(), flip.end(),
    //   boost::variate_generator<rng_type&, boost::uniform_smallint<> >(
    //   uniform_01, boost::uniform_smallint<>(0, 1)));

    if (param.model.is_signed()) {
      int n = 0;
      std::vector<int>::iterator itr = flip.begin();
      std::vector<int>::iterator itr_end = flip.end();
      std::vector<int>::iterator sitr = config.loop_sign.begin();
      for (; itr != itr_end; ++itr, ++sitr) if (*itr == 1 && *sitr == -1) ++n;
      if (n % 2 == 1) config.sign = -config.sign;
    }

    // flip spins
    {
      vertex_iterator vi, vi_end;
      for (boost::tie(vi, vi_end) = boost::vertices(param.vgraph);
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
      for (boost::tie(vi, vi_end) = boost::vertices(param.vgraph);
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

  // measurements

  static int loop_index_0(int i, const config_type& config)
  { return config.wl.series(i).first->loop_segment(0).index; }

  static int loop_index_1(int i, const config_type& config)
  { return config.wl.series(i).second->loop_segment(0).index; }

  static double static_sz(int i, const config_type& config)
  { return 0.5 - double(config.wl.series(i).first->conf()); }

  static double dynamic_sz(int i, const config_type& config, const graph_type&)
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
