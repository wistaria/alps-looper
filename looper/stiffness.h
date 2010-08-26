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

#ifndef LOOPER_STIFFNESS_H
#define LOOPER_STIFFNESS_H

#include "measurement.h"

namespace looper {

//
// stiffness measurement
//

template<unsigned int MAX_DIM>
struct stiffness : public has_improved_estimator_tag, public has_normal_estimator_tag {

  template<typename MC, typename LAT, typename TIME>
  struct estimator
  {
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<real_bond_t,
      const typename lattice_t::virtual_graph_type,
      typename real_bond_descriptor<lattice_t>::type>::type
      real_bond_map_t;
    typedef typename alps::property_map<alps::bond_vector_relative_t,
      const typename lattice_t::real_graph_type,
      alps::coordinate_type>::type bond_vector_relative_map_t;

    real_bond_map_t real_bond;
    bond_vector_relative_map_t bond_vector_relative;
    unsigned int dim;

    void initialize(alps::Parameters const& /* params */, lattice_t const& lat,
      bool /* is_signed */, bool enable_improved_estimator) {
      real_bond = alps::get_or_default(real_bond_t(), lat.vg(),
        typename real_bond_descriptor<lattice_t>::type());
      bond_vector_relative = alps::get_or_default(bond_vector_relative_t(), lat.rg(),
        coordinate_type());
      dim = get_property(lat.rg(), dimension_t());
      if (enable_improved_estimator && dim > MAX_DIM) {
        std::cerr << "Spatial dimension (=" << dim << ") is too large.  "
                  << "Stiffness will be measured only for the first "
                  << MAX_DIM << " dimensions\n";
        dim = MAX_DIM;
      }
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool) {
      if (dim > 0) add_scalar_obs(m, "Stiffness", is_signed);
    }

    struct improved_estimator {
      struct estimate {
        double winding[MAX_DIM];
        estimate() : winding() {}
        void reset(estimator_t const&) {
          for (int i = 0; i < MAX_DIM; ++i) winding[i] = 0;
        }
        estimate& operator+=(estimate const& rhs) {
          for (int i = 0; i < MAX_DIM; ++i) winding[i] += rhs.winding[i];
          return *this;
        }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) const {}
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double, int b, int, int c) {
          alps::coordinate_type const& vr =
            emt.bond_vector_relative[emt.real_bond[bond(lat.vg(), b)]];
          for (int i = 0; i < emt.dim; ++i) winding[i] -= (1-2*c) * vr[i];
        }
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) const {}
        void end_bs(estimator_t const& emt, lattice_t const& lat, double, int b, int, int c) {
          alps::coordinate_type const& vr =
            emt.bond_vector_relative[emt.real_bond[bond(lat.vg(), b)]];
          for (int i = 0; i < emt.dim; ++i) winding[i] += (1-2*c) * vr[i];
        }
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const&, lattice_t const&, double, int, int) const {}
        void start(estimator_t const&, lattice_t const&, double, int, int) const {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) const {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) const {}
      };

      struct collector {
        unsigned int dim;
        double w2;
        void reset(estimator_t const& emt) { dim = emt.dim; w2 = 0; }
        collector& operator+=(collector const& rhs) {
          for (int i = 0; i < dim; ++i) w2 += rhs.w2;
          return *this;
        }
        collector& operator+=(estimate const& rhs) {
          for (int i = 0; i < dim; ++i) w2 += power2(0.5 * rhs.winding[i]);
          return *this;
        }
        template<typename M>
        void commit(M& m, lattice_t const&, double beta, double sign, double) const {
          if (dim > 0) m["Stiffness"] << sign * w2 / (beta * dim);
        }
      };
    };

    struct normal_estimator {
      struct collector {
        unsigned int dim;
        double winding[MAX_DIM];
        void reset(estimator_t const& emt) {
          dim = emt.dim;
          for (int i = 0; i < MAX_DIM; ++i) winding[i] = 0;
        }
        collector& operator+=(collector const& rhs) {
          for (int i = 0; i < MAX_DIM; ++i) winding[i] += rhs.winding[i];
          return *this;
        }
        void begin_s(estimator_t const&, lattice_t const&, double, int, int) const {}
        void begin_bs(estimator_t const& emt, lattice_t const& lat, double, int b, int, int c) {
          alps::coordinate_type const& vr =
            emt.bond_vector_relative[emt.real_bond[bond(lat.vg(), b)]];
          for (int i = 0; i < emt.dim; ++i) winding[i] -= (1-2*c) * vr[i];
        }
        void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void end_s(estimator_t const&, lattice_t const&, double, int, int) const {}
        void end_bs(estimator_t const& emt, lattice_t const& lat, double, int b, int, int c) {
          alps::coordinate_type const& vr =
            emt.bond_vector_relative[emt.real_bond[bond(lat.vg(), b)]];
          for (int i = 0; i < emt.dim; ++i) winding[i] += (1-2*c) * vr[i];
        }
        void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
        void start_bottom(estimator_t const&, lattice_t const&, double, int, int) const {}
        void start(estimator_t const&, lattice_t const&, double, int, int) const {}
        void stop(estimator_t const&, lattice_t const&, double, int, int) const {}
        void stop_top(estimator_t const&, lattice_t const&, double, int, int) const {}
        template<typename M>
        void commit(M& m, lattice_t const&, double beta, double sign, double) const {
          double w2 = 0;
          for (int i = 0; i < dim; ++i) w2 += power2(0.5 * winding[i]);
          if (dim > 0) m["Stiffness"] << sign * w2 / (beta * dim);
        }
      };
    };
  };
};

} // end namespace looper

#endif // LOOPER_STIFFNESS_H
