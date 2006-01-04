/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_CLUSTER_H
#define LOOPER_CLUSTER_H

#include "lattice.h"
#include "type.h"
#include <boost/mpl/bool.hpp>

namespace looper {

struct cluster_info
{
  cluster_info(bool t = false) : to_flip(t), weight(0) {}
  bool to_flip;
  double weight;
};

struct cluster_measure
{
  double usize0;
  double ssize0;
  double umag0;
  double smag0;
  double usize;
  double ssize;
  double umag;
  double smag;
  cluster_measure()
    : usize0(0), ssize0(0), umag0(0), smag0(0),
      usize(0), ssize(0), umag(0), smag(0)
  {}

  template<typename G>
  void at_zero(boost::mpl::true_ const&, G const& g, int s, double c)
  {
    usize0 += 0.5;
    ssize0 += 0.5 * gauge(g, s);
    umag0  += c;
    smag0  += c * gauge(g, s);
  }
  template<typename G>
  void at_zero(boost::mpl::false_ const&, G const&, int, double c)
  {
    usize0 += 0.5;
    umag0  += c;
  }

  template<typename G>
  void start(boost::mpl::true_ const&, G const& g, int s, double t, double c)
  {
    usize -= t;
    ssize -= t * gauge(g, s);
    umag  -= c * t;
    smag  -= c * t * gauge(g, s);
  }
  template<typename G>
  void start(boost::mpl::false_ const&, G const&, int, double t, double c)
  {
    usize -= t;
    umag  -= c * t;
  }

  template<typename G>
  void term(boost::mpl::true_ const&, G const& g, int s, double t, double c)
  {
    usize += t;
    ssize += t * gauge(g, s);
    umag  += c * t;
    smag  += c * t * gauge(g, s);
  }
  template<typename G>
  void term(boost::mpl::false_ const&, G const&, int, double t, double c)
  {
    usize += t;
    umag  += c * t;
  }
};

struct cluster_accumulator
{
  double usize0;
  double ssize0;
  double umag0;
  double smag0;
  double usize;
  double ssize;
  double umag;
  double smag;
  cluster_accumulator()
    : usize0(0), ssize0(0), umag0(0), smag0(0),
      usize(0), ssize(0), umag(0), smag(0)
  {}

  cluster_accumulator operator+(cluster_measure const& cm)
  {
    using looper::sqr;
    usize0 += sqr(cm.usize0);
    ssize0 += sqr(cm.ssize0);
    umag0  += sqr(cm.umag0);
    smag0  += sqr(cm.smag0);
    usize  += sqr(cm.usize);
    ssize  += sqr(cm.ssize);
    umag   += sqr(cm.umag);
    smag   += sqr(cm.smag);
    return *this;
  }

  template<class M>
  void commit(M& measurements, bool bipartite, int n) const
  {
    measurements["Magnetization"] << 0.0;
    measurements["Magnetization Density"] << 0.0;
    measurements["Magnetization^2"] << umag0;
    measurements["Susceptibility"] << umag / n;
    measurements["Generalized Magnetization^2"] << usize0;
    measurements["Generalized Susceptibility"] << usize / n;
    if (bipartite) {
      measurements["Staggered Magnetization"] << 0.0;
      measurements["Staggered Magnetization Density"] << 0.0;
      measurements["Staggered Magnetization^2"] << smag0;
      measurements["Staggered Susceptibility"] << smag / n;
      measurements["Staggered Generalized Magnetization^2"] << ssize0;
      measurements["Staggered Generalized Susceptibility"] << ssize / n;
    }
  }
};

} // end namespace looper

#endif // LOOPER_CLUSTER_H
