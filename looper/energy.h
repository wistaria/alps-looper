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

#ifndef LOOPER_ENERGY_H
#define LOOPER_ENERGY_H

#include "measurement.h"

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct energy : public has_normal_estimator_tag, public has_evaluator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;

    void initialize(alps::Parameters const&, lattice_t const&, bool /* is_signed */,
      bool /* enable_improved_estimator */) {}
    template<typename M>
    void init_observables(M& m, bool is_signed, bool /* enable_improved_estimator */) {
      add_scalar_obs(m, "Energy", is_signed);
      add_scalar_obs(m, "Energy Density", is_signed);
      add_scalar_obs(m, "Energy^2", is_signed);
    }

    struct normal_estimator {
      struct collector {
        double ene_;
        collector() : ene_(0) {}
        void reset(estimator_t const&) { ene_ = 0; }
        collector& operator+=(collector const& rhs) {
          ene_ += rhs.ene_;
          return *this;
        }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
        void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const&, lattice_t const&, double, int, int) {}
        void start(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) {}
        void set_energy(double ene) { ene_ = ene; }
        template<typename M>
        void commit(M& m, estimator_t const&, lattice_t const& lat, double beta, double sign,
          double nop, std::vector<int> const&) const {
          m["Energy"] << sign * ene_;
          m["Energy Density"] << sign * ene_ / lat.volume();
          m["Energy^2"] << sign * (power2(ene_) - nop / power2(beta));
#ifdef STD_OUTPUT
          std::cout << sign * ene_  / lat.volume() << ' ';
#endif
        }
      };
    };
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const&,
      alps::ObservableSet const& m_in) {
      try {
        alps::RealObsevaluator beta = m_in["Inverse Temperature"];
        alps::RealObsevaluator vol = m_in["Volume"];
        alps::RealObsevaluator ene = m_in["Energy"];
        alps::RealObsevaluator ene2 = m_in["Energy^2"];
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / vol.mean();
        m.addObservable(c);
      } catch (...) {}
    }
  };
};

} // end namespace looper

#endif // LOOPER_ENERGY_H
