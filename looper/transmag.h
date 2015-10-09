/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006-2011 by Synge Todo <wistaria@comp-phys.org>
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

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct transverse_magnetization : public has_improved_estimator_tag {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;

    void initialize(alps::Parameters const& /* params */, lattice_t const& /* lat */,
      bool /* is_signed */, bool enable_improved_estimator) {
      if (!enable_improved_estimator)
        std::cout << "Warning: transverse magnetization measurement is disabled.\n";
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      if (enable_improved_estimator) {
        add_scalar_obs(m, "Transverse Magnetization", is_signed);
        add_scalar_obs(m, "Transverse Magnetization Density", is_signed);
      }
    }

    struct improved_estimator {
      struct estimate {
        double length;
        bool closed;
        estimate() : length(0), closed(true) {}
        void reset(estimator_t const&) {
          length = 0;
          closed = true;
        }
        estimate& operator+=(estimate const& rhs) {
          length += rhs.length;
          closed &= rhs.closed;
          return *this;
        }
        void begin_s(estimator_t const&, lattice_t const&, double t, int, int) {
          length -= t;
          closed = false;
        }
        void begin_bs(estimator_t const&, lattice_t const&, double t, int, int, int) {
          length -= t;
        }
        void begin_bt(estimator_t const&, lattice_t const&, double t, int, int, int) {
          length -= t;
        }
        void end_s(estimator_t const&, lattice_t const&, double t, int, int) {
          length += t;
          closed = false;
        }
        void end_bs(estimator_t const&, lattice_t const&, double t, int, int, int) { length += t; }
        void end_bt(estimator_t const&, lattice_t const&, double t, int, int, int) { length += t; }
        void start_bottom(estimator_t const&, lattice_t const&, double t, int, int) {
          length -= t;
        }
        void start(estimator_t const&, lattice_t const&, double t, int, int) { length -= t; }
        void stop(estimator_t const&, lattice_t const&, double t, int, int) { length += t; }
        void stop_top(estimator_t const&, lattice_t const&, double t, int, int) { length += t; }
      };

      struct collector {
        double length;
        void reset(estimator_t const&) { length = 0; }
        collector& operator+=(collector const& coll) { length += coll.length; }
        collector& operator+=(estimate const& est) {
          length += (est.closed ? 0.0 : est.length);
          return *this;
        }
        template<typename M>
        void commit(M& m, estimator_t const&, lattice_t const& lat, double, double sign, double,
          std::vector<int> const&) const {
          m["Transverse Magnetization"] << 0.5 * sign * length;
          m["Transverse Magnetization Density"] << 0.5 * sign * length / lat.volume();
        }
      };
    };
  };
};

} // end namespace looper

#endif // TRANSMAG_H
