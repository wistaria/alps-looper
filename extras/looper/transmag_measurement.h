/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef TRANSMAG_MEASUREMENT_H
#define TRANSMAG_MEASUREMENT_H

#include <looper/measurement.h>
#include <alps/alea.h>
#include <cmath>
#include <string>

struct transverse_magnetization_estimator : public looper::base_estimator
{
  template<typename T>
  static void initialize(T& m, bool, bool is_signed,
                         bool use_improved_estimator)
  {
    if (use_improved_estimator) {
      looper::add_measurement(m, "Transverse Magnetization",
                              is_signed);
      looper::add_measurement(m, "Transverse Magnetization Density",
                              is_signed);
    }
  }

  struct estimate : public looper::base_estimator::estimate
  {
    double length;
    bool closed;
    estimate() : length(0), closed(true) {}
    template<typename G>
    void start1(G const&, double t, int, int)
    {
      length -= t;
      closed = false;
    }
    template<typename G>
    void start2(G const&, double t, int, int)
    { length -= t; }
    template<typename G>
    void term1(G const&, double t, int, int)
    {
      length += t;
      closed = false;
    }
    template<typename G>
    void term2(G const&, double t, int, int)
    { length += t; }
    template<typename G>
    void at_bot(G const&, double t, int, int)
    { length -= t; }
    template<typename G>
    void at_top(G const&, double t, int, int)
    { length += t; }
  };

  template<typename QMC, typename IMPROVE>
  struct collector
  {
    double length;
    collector() : length(0) {}
    template<typename EST>
    collector operator+(EST const& cm)
    {
      if (!cm.closed) length += cm.length;
      return *this;
    }
    template<typename M>
    void commit(M& m, bool, double, int nrs, int, double sign) const
    {
      m["Transverse Magnetization"] << 0.5 * sign * length;
      m["Transverse Magnetization Density"] << 0.5 * sign * length / nrs;
    }
  };
};

#endif // TRANSMAG_MEASUREMENT_H
