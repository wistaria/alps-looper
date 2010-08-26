/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_DUMB_MEASUREMENT_H
#define LOOPER_DUMB_MEASUREMENT_H

#include "type.h"
#include <alps/parameter/parameters.h>
#include <alps/alea/observableset.h>
#include <vector>

namespace looper {

template<typename MEASUREMENT>
struct dumb_measurement : public has_improved_estimator_tag, public has_normal_estimator_tag,
  public has_evaluator_tag {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef estimator<MC, LAT, TIME> estimator_t;
    typedef LAT lattice_t;
    void initialize(alps::Parameters const&, lattice_t const& lat, bool /* is_signed */) {}
    template<typename M>
    void init_observables(M&, bool /* is_signed */, bool /* enable_improved_estimator */) {}

    struct improved_estimator {
      struct estimate {
        void reset() {}
        estimate& operator+=(estimate const&) { return *this; }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int) {}
        void start_bottom(estimator_t const&, lattice_t const&, double, int, int) {}
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
      };
      struct collector {
        collector() {}
        void reset() {}
        collector& operator+=(collector const&) { return *this; }
        collector& operator+=(estimate const&) { return *this; }
        template<typename M>
        void commit(M&, lattice_t const&, double /* beta */, double /* sign */,
          int /* nop */) const {}
      };
    };

    struct normal_estimator {
      struct collector {
        collector() {}
        void reset() {}
        collector& operator+=(collector const&) { return *this; }
        void begin_s(estimator_t const&, lattice_t const&, double /* t */, int /* s */,
          int /* c */) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int) {}
        void start_bottom(estimator_t const&, lattice_t const&, double, int, int) {}
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
        template<typename M>
        void commit(M&, lattice_t const&, double /* beta */, double /* sign */,
          int /* nop */) const {}
      };
    };
  };

  struct pre_evaluator {
    static void pre_evaluate(alps::ObservableSet&, alps::Parameters const&,
      alps::ObservableSet const&) {}
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet&, alps::Parameters const&,
      alps::ObservableSet const&) {}
  };
};

} // end namespace looper

#endif // LOOPER_DUMB_MEASUREMENT_H
