/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: measurement.h 554 2003-11-12 02:36:24Z wistaria $

#ifndef LOOPER_MEASUREMENT_H
#define LOOPER_MEASUREMENT_H

namespace looper {

template<class T> T sqr(T t) { return t * t; }

//
// unimproved estimators
//

// energy offset

template<class P>
inline double energy_offset(const P& param)
{
  return P::qmc_type::energy_offset(param.virtual_graph, param.model);
}

template<class C, class P>
inline double energy_z(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;

  double ene = 0.;
  typename qmc_type::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(param.virtual_graph.graph);
       ei != ei_end; ++ei)
    ene -= param.model.bond(bond_type(*ei, param.virtual_graph.graph)).jz() *
      qmc_type::static_sz(boost::source(*ei, param.virtual_graph.graph),
			  config) *
      qmc_type::static_sz(boost::target(*ei, param.virtual_graph.graph),
			  config);
  return ene / param.virtual_graph.num_real_vertices;
}


// energy

template<class C, class P>
inline double energy_xy(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;
  return - (double)qmc_type::num_offdiagonals(config) / param.beta /
    param.virtual_graph.num_real_vertices;
}

template<class C, class P>
inline boost::tuple<double, double, double>
energy(const C& config, const P& param)
{
  typedef typename C::qmc_type qmc_type;
  double ez = energy_z(config, param);
  double exy = energy_xy(config, param);
  double e2 = sqr(ez + exy) -
    (double)qmc_type::num_offdiagonals(config) /
    sqr(param.beta * param.virtual_graph.num_real_vertices);
  return boost::make_tuple(ez, exy, e2);
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
    double sd = 0.;
    double ss = 0.;
    typename qmc_type::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
	 vi != vi_end; ++vi) {
      sd += gauge(*vi, param.virtual_graph.graph) *
	qmc_type::dynamic_sz(*vi, config, param.virtual_graph);
      ss += gauge(*vi, param.virtual_graph.graph) *
	qmc_type::static_sz(*vi, config);
    }
    //// return param.beta * (sqr(sd) + sqr(ss) / (config.os.size() + 1)) /
    //// param.virtual_graph.num_real_vertices;
    return param.beta * sqr(sd) /
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
    std::vector<double> len(config.num_loops);
    std::vector<double> mg(config.num_loops0);
    typename qmc_type::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(param.virtual_graph.graph);
	 vi != vi_end; ++vi) {
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

      ++mg[qmc_type::loop_index_0(*vi, config)];
    }      

    double ss = 0.;
    for (int i = 0; i < config.num_loops; ++i) ss += sqr(len[i]);
    ss *= 0.25;
    double m2 = 0.;
    for (int i = 0; i < config.num_loops0; ++i) m2 += sqr(mg[i]);
    m2 *= 0.25;

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
    std::vector<double> len(config.num_loops);
    std::vector<double> mg(config.num_loops0);

    std::vector<int> pos(boost::num_vertices(param.virtual_graph.graph), 0);
    int p = 1;
    typename qmc_type::config_type::const_iterator oi_end = config.os.end();
    for (typename qmc_type::config_type::const_iterator oi = config.os.begin();
	 oi != oi_end; ++oi, ++p) {
      if (oi->is_offdiagonal()) {
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

    double ss = 0.;
    for (int i = 0; i < config.num_loops; ++i) ss += sqr(len[i]);
    ss *= 0.25 / sqr(config.os.size());
    double m2 = 0.;
    for (int i = 0; i < config.num_loops0; ++i) m2 += sqr(mg[i]);
    m2 *= 0.25;

    return boost::make_tuple(m2 / param.virtual_graph.num_real_vertices,
			     param.beta *
			     (ss /* +
			       m2 / (config.os.size() + 1) */) /
			     param.virtual_graph.num_real_vertices);
//     return boost::make_tuple(0.25 * m2 / param.virtual_graph.num_real_vertices,
// 			     0.25 * param.beta *
// 			     (ss / config.os.size() +
// 			      m2 / (config.os.size() + 1)) /
// 			     param.virtual_graph.num_real_vertices);
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
