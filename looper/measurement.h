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

#ifndef LOOPER_MEASUREMENT_H
#define LOOPER_MEASUREMENT_H

#include <looper/util.h>

namespace looper {

//
// unimproved estimators
//

// energy offset

template<class P>
inline double energy_offset(const P& param)
{
  return P::qmc_type::energy_offset(param.virtual_graph, param.model);
}


// energy

namespace {

template<class Q> struct energy_helper;

template<class G, class M, class W, class N>
struct energy_helper<path_integral<G, M, W, N> >
{
  typedef path_integral<G, M, W, N> qmc_type;
  static boost::tuple<double, double, double>
  calc(const typename qmc_type::config_type& config,
       const typename qmc_type::parameter_type& param)
  {
    double ez = 0.;
    typename qmc_type::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(param.virtual_graph.graph);
         ei != ei_end; ++ei) {
      ez -= param.model.bond(bond_type(*ei, param.virtual_graph.graph)).jz()
        * qmc_type::static_correlation(config,
           boost::source(*ei, param.virtual_graph.graph),
           boost::target(*ei, param.virtual_graph.graph));
    }
    double exy = -(double)config.wl.num_links() / param.beta;
    double e2 = sqr(ez + exy) -
      (double)config.wl.num_links() / sqr(param.beta);
    double nvinv = 1.0 / (double)param.virtual_graph.num_real_vertices;
    return boost::make_tuple(nvinv * ez, nvinv * exy, sqr(nvinv) * e2);
  }
};

template<class G, class M, class W, class N>
struct energy_helper<sse<G, M, W, N> >
{
  typedef sse<G, M, W, N> qmc_type;
  static boost::tuple<double, double, double>
  calc(const typename qmc_type::config_type& config,
       const typename qmc_type::parameter_type& param)
  {
    typedef typename qmc_type::config_type::const_iterator
      const_operator_iterator;
    int nd = 0;
    int no = 0;
    const_operator_iterator oi_end = config.os.end();
    for (const_operator_iterator oi = config.os.begin(); oi != oi_end; ++oi) {
      if (oi->is_diagonal()) ++nd;
      if (oi->is_offdiagonal()) ++no;
    }
    double ez = - double(nd) / param.beta - param.ez_offset;
    double exy = - double(no) / param.beta;
    double e2 = sqr(ez + exy) - (double)(nd + no) / sqr(param.beta);
    double nvinv = 1.0 / (double)param.virtual_graph.num_real_vertices;
    return boost::make_tuple(nvinv * ez, nvinv * exy, sqr(nvinv) * e2);
  }
};

} // end namespace

template<class C, class P>
inline boost::tuple<double, double, double>
energy(const C& config, const P& param)
{
  return energy_helper<typename C::qmc_type>::calc(config, param);
}


// total Sz

template<class C, class P>
inline double uniform_sz(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;

  double sz = 0.;
  typename qmc_type::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
       vi != vi_end; ++vi) {
    sz += qmc_type::static_sz(*vi, config);
  }
  return sz / param.virtual_graph.num_real_vertices;
}


// total stsaggered Sz

template<class C, class P>
inline double staggered_sz(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;
  if (!param.is_bipartite)
    boost::throw_exception(std::runtime_error("lattice is not bipartitle"));
  double ss = 0.;
  typename qmc_type::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
       vi != vi_end; ++vi)
    ss += gauge(*vi, param.virtual_graph.graph) *
      qmc_type::static_sz(*vi, config);
  return ss / param.virtual_graph.num_real_vertices;
}


// staggered susceptilibity

namespace {

template<class Q> struct staggered_susceptibility_helper;

template<class G, class M, class W, class N>
struct staggered_susceptibility_helper<path_integral<G, M, W, N> >
{
  typedef path_integral<G, M, W, N> qmc_type;
  static double calc(const typename qmc_type::config_type& config,
                     const typename qmc_type::parameter_type& param)
  {
    double ss = 0.;
    typename qmc_type::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
         vi != vi_end; ++vi)
      ss += gauge(*vi, param.virtual_graph.graph) *
        qmc_type::dynamic_sz(*vi, config, param.virtual_graph);
    return param.beta * ss * ss / param.virtual_graph.num_real_vertices;
  }
};

template<class G, class M, class W, class N>
struct staggered_susceptibility_helper<sse<G, M, W, N> >
{
  typedef sse<G, M, W, N> qmc_type;
  static double calc(const typename qmc_type::config_type& config,
                     const typename qmc_type::parameter_type& param)
  {
    std::vector<double> conf(boost::num_vertices(param.virtual_graph.graph));
    double ss = 0.;
    typename qmc_type::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
         vi != vi_end; ++vi) {
      double c = gauge(*vi, param.virtual_graph.graph) *
	qmc_type::static_sz(*vi, config);
      conf[*vi] = c;
      ss += c;
    }
    double sst = ss;
    double sd = 0.;
    typename qmc_type::config_type::const_iterator oi_end = config.os.end();
    for (typename qmc_type::config_type::const_iterator oi = config.os.begin();
	 oi != oi_end; ++oi) {
      sd += sst;
      if (oi->is_offdiagonal()) {
	typename qmc_type::edge_iterator ei =
	  boost::edges(param.virtual_graph.graph).first + oi->bond();
	int s0 = boost::source(*ei, param.virtual_graph.graph);
	conf[s0] *= -1;
	sst += (2 * conf[s0]);
	int s1 = boost::target(*ei, param.virtual_graph.graph);
	conf[s1] *= -1;
	sst += (2 * conf[s1]);
      }
    }
    sd /= std::sqrt((double)config.os.size() * (double)(config.os.size() + 1));
    return param.beta * (sqr(sd) + sqr(ss) / (config.os.size() + 1)) /
      param.virtual_graph.num_real_vertices;
  }
};

} // end namespace

template<class C, class P>
inline double staggered_susceptibility(const C& config, const P& param)
{
  if (!param.is_bipartite)
    boost::throw_exception(std::runtime_error("lattice is not bipartitle"));
  return staggered_susceptibility_helper<typename C::qmc_type>::
    calc(config, param);
}


//
// improved estimators
//

// diagonal energy

template<class C, class P>
inline double energy_z_imp(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;
  double ene = 0.;
  typename qmc_type::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(param.virtual_graph.graph);
       ei != ei_end; ++ei) {
    typename qmc_type::vertex_descriptor v0 =
      boost::source(*ei, param.virtual_graph.graph);
    typename qmc_type::vertex_descriptor v1 =
      boost::target(*ei, param.virtual_graph.graph);
    if (qmc_type::loop_index_0(v0, config) ==
        qmc_type::loop_index_0(v1, config))
      ene -= param.model.bond(bond_type(*ei, param.virtual_graph.graph)).jz() *
        qmc_type::static_sz(v0, config) * qmc_type::static_sz(v1, config);
  }
  return ene / param.virtual_graph.num_real_vertices;
}

// uniform magnetizations^2

template<class C, class P>
inline double
uniform_sz2_imp(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;

  std::vector<double> m(config.num_loops0);
  typename qmc_type::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
       vi != vi_end; ++vi)
    m[qmc_type::loop_index_0(*vi, config)] += qmc_type::static_sz(*vi, config);

  double m2 = 0.;
  for (int i = 0; i < config.num_loops0; ++i) m2 += sqr(m[i]);

  return m2 / param.virtual_graph.num_real_vertices;
}

// staggered magnetizations^2

template<class C, class P>
inline double
staggered_sz2_imp(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;

  std::vector<double> m(config.num_loops0);
  typename qmc_type::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
       vi != vi_end; ++vi)
    m[qmc_type::loop_index_0(*vi, config)] +=
      gauge(*vi, param.virtual_graph.graph) * qmc_type::static_sz(*vi, config);

  double m2 = 0.;
  for (int i = 0; i < config.num_loops0; ++i) m2 += sqr(m[i]);

  return m2 / param.virtual_graph.num_real_vertices;
}


// generalized susceptilibity

namespace {

template<class Q> struct generalized_susceptibility_imp_helper;

template<class G, class M, class W, class N>
struct generalized_susceptibility_imp_helper<path_integral<G, M, W, N> >
{
  typedef path_integral<G, M, W, N> qmc_type;
  static boost::tuple<double, double>
  calc(const typename qmc_type::config_type& config,
       const typename qmc_type::parameter_type& param)
  {
    std::vector<double> mg(config.num_loops0);
    std::vector<double> len(config.num_loops);
    typename qmc_type::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
         vi != vi_end; ++vi) {
      ++mg[qmc_type::loop_index_0(*vi, config)];

      // setup iterators
      typename qmc_type::config_type::const_iterator
        itrD = config.wl.series(*vi).first;
      typename qmc_type::config_type::const_iterator
        itrU = boost::next(itrD);

      // iteration up to t = 1
      for (;; itrD = itrU, ++itrU) {
        len[qmc_type::segment_d(itrU).index] += (itrU->time() - itrD->time());
        if (itrU.at_top()) break;
      }
    }

    double m2 = 0.;
    for (int i = 0; i < config.num_loops0; ++i) m2 += sqr(mg[i]);
    m2 *= 0.25;
    double ss = 0.;
    for (int i = 0; i < config.num_loops; ++i) ss += sqr(len[i]);
    ss *= 0.25;

    return boost::make_tuple(m2 / param.virtual_graph.num_real_vertices,
                             param.beta * ss /
                             param.virtual_graph.num_real_vertices);
  }
};

template<class G, class M, class W, class N>
struct generalized_susceptibility_imp_helper<sse<G, M, W, N> >
{
  typedef sse<G, M, W, N> qmc_type;
  static boost::tuple<double, double>
  calc(const typename qmc_type::config_type& config,
       const typename qmc_type::parameter_type& param)
  {
    std::vector<double> mg(config.num_loops0);
    std::vector<double> len(config.num_loops);

    std::vector<int> pos(boost::num_vertices(param.virtual_graph.graph), 0);
    int p = 1;
    typename qmc_type::config_type::const_iterator oi_end = config.os.end();
    for (typename qmc_type::config_type::const_iterator oi = config.os.begin();
         oi != oi_end; ++oi, ++p) {
      if (!oi->is_identity()) {
        typename qmc_type::edge_iterator ei =
          boost::edges(param.virtual_graph.graph).first + oi->bond();
        typename qmc_type::vertex_descriptor v0 =
          boost::source(*ei, param.virtual_graph.graph);
        typename qmc_type::vertex_descriptor v1 =
          boost::target(*ei, param.virtual_graph.graph);
        len[qmc_type::segment_d(oi, 0).index] += (p - pos[v0]);
        len[qmc_type::segment_d(oi, 1).index] += (p - pos[v1]);
        pos[v0] = p;
        pos[v1] = p;
      }
    }
    typename qmc_type::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
         vi != vi_end; ++vi) {
      len[qmc_type::loop_index_1(*vi, config)] +=
        (config.os.size() - pos[*vi]);
      ++mg[qmc_type::loop_index_0(*vi, config)];
    }

    double m2 = 0.;
    for (int i = 0; i < config.num_loops0; ++i) m2 += sqr(mg[i]);
    m2 *= 0.25;
    double ss = 0.;
    for (int i = 0; i < config.num_loops; ++i) ss += sqr(len[i]);
    ss *= 0.25 / ((double)config.os.size() * (double)(config.os.size() + 1));

    return boost::make_tuple(m2 / param.virtual_graph.num_real_vertices,
      param.beta * (ss + m2 / (double)(config.os.size() + 1)) /
      param.virtual_graph.num_real_vertices);
  }
};

} // end namespace

template<class C, class P>
inline boost::tuple<double, double>
generalized_susceptibility_imp(const C& config, const P& param)
{
  return generalized_susceptibility_imp_helper<typename C::qmc_type>::
    calc(config, param);
}

} // namespace looper

#endif // LOOPER_MEASUREMENT_H
