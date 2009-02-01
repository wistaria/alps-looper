/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TRANSMAG_H
#define LOOPER_TRANSMAG_H

#ifndef LOOPER_ONLY_PATH_INTEGRAL
# define LOOPER_ONLY_PATH_INTEGRAL
#endif

#include "measurement.h"

namespace looper {

struct transverse_magnetization {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;

    template<typename M>
    void initialize(M& m, alps::Parameters const& /* params */, lattice_t const& /* lat */,
      bool is_signed, bool use_improved_estimator) {
      if (use_improved_estimator) {
        add_scalar_obs(m, "Transverse Magnetization", is_signed);
        add_scalar_obs(m, "Transverse Magnetization Density", is_signed);
      }
    }

    // improved estimator

    struct estimate {
      double length;
      bool closed;
      estimate() : length(0), closed(true); {}
      void init() {
        length = 0;
        closed = true;
      }
      void begin_s(lattice_t const&, double t, int, int) {
        length -= t;
        closed = false;
      }
      void begin_bs(lattice_t const&, double t, int, int, int) { length -= t; }
      void begin_bt(lattice_t const&, double t, int, int, int) { length -= t; }
      void end_s(lattice_t const&, double t, int, int) {
        length += t;
        closed = false;
      }
      void end_bs(lattice_t const&, double t, int, int, int) { length += t; }
      void end_bt(lattice_t const&, double t, int, int, int) { length += t; }
      void start_bottom(lattice_t const&, double t, int, int) { length -= t; }
      void start(lattice_t const&, double t, int, int) { length -= t; }
      void stop(lattice_t const&, double t, int, int) { length += t; }
      void stop_top(lattice_t const&, double t, int, int) { length += t; }
    };
    void init_estimate(estimate& es) const { es.init(); }

    struct collector {
      double length;
      void init() { length = 0; }
      collector& operator+=(estimate const& cm) {
        if (!cm.closed) length += cm.length;
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const& lat, double, int, double sign) const {
        m["Transverse Magnetization"] << 0.5 * sign * length;
        m["Transverse Magnetization Density"] << 0.5 * sign * length / lat.volume();
      }
    };
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& /* spins */, std::vector<OP> const& operators,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& /* m */, lattice_t const& /* lat */, double /* beta */,
      double /* sign */, std::vector<int> const& /* spins */,
      std::vector<OP> const& /* operators */, std::vector<int> const& /* spins_c */) {}
  };
};

} // end namespace looper

#endif // TRANSMAG_H
