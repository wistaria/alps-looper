/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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
#include <alps/fixed_capacity_vector.h>

namespace looper {

//
// stiffness measurement
//

template<unsigned int MAX_DIM>
struct stiffness
{
  template<typename MC, typename VLAT, typename TIME>
  struct estimator
  {
    typedef MC   mc_type;
    typedef VLAT virtual_lattice_t;
    typedef TIME time_t;
    typedef typename alps::property_map<alps::bond_vector_relative_t,
              const typename virtual_lattice_t::real_graph_type,
              alps::coordinate_type>::type bond_vector_relative_map_t;

    bool improved;
    bond_vector_relative_map_t bond_vector_relative;
    unsigned int dim;

    template<typename M>
    void initialize(M& m, alps::Parameters const& /* params */,
                    virtual_lattice_t const& vlat,
                    bool is_signed, bool use_improved_estimator)
    {
      improved = use_improved_estimator;
      bond_vector_relative =
        get_or_default(alps::bond_vector_relative_t(), vlat.rgraph(),
                       alps::coordinate_type());
      dim = boost::get_property(vlat.rgraph(), dimension_t());
      if (improved && dim > MAX_DIM) {
        std::cerr << "Spatial dimension (=" << dim << ") is too large.  "
                  << "Stiffness will be measured only for the first "
                  << MAX_DIM << " dimensions\n";
        dim = MAX_DIM;
      }
      if (dim > 0) add_scalar_obs(m, "Stiffness", is_signed);
    }

    // improved estimator

    struct estimate
    {
      bond_vector_relative_map_t bond_vector_relative;
      alps::fixed_capacity_vector<double, MAX_DIM> winding;
      void init(bond_vector_relative_map_t map, int dim)
      {
        bond_vector_relative = map;
        winding.resize(dim);
      }
      void start_s(virtual_lattice_t const&, double, int, int) const {}
      void start_bs(virtual_lattice_t const& vlat, double, int b, int, int c)
      {
        alps::coordinate_type const& vr = bond_vector_relative[rbond(vlat, b)];
        for (int i = 0; i < winding.size(); ++i) winding[i] -= (1-2*c) * vr[i];
      }
      void start_bt(virtual_lattice_t const&, double, int, int, int) {}
      void term_s(virtual_lattice_t const&, double, int, int) const {}
      void term_bs(virtual_lattice_t const& vlat, double, int b, int,
                   int c)
      {
        alps::coordinate_type const& vr = bond_vector_relative[rbond(vlat, b)];
        for (int i = 0; i < winding.size(); ++i) winding[i] += (1-2*c) * vr[i];
      }
      void term_bt(virtual_lattice_t const&, double, int, int, int) {}
      void at_bot(virtual_lattice_t const&, double, int, int) const {}
      void at_top(virtual_lattice_t const&, double, int, int) const {}
    };
    void init_estimate(estimate& est) const
    { est.init(bond_vector_relative, dim); }

    struct collector
    {
      unsigned int dim;
      double w2;
      void init(unsigned int d) { dim = d; w2 = 0; }
      template<typename EST>
      collector operator+(EST const& est)
      {
        for (int i = 0; i < dim; ++i) w2 += power2(0.5 * est.winding[i]);
        return *this;
      }
      template<typename M>
      void commit(M& m, virtual_lattice_t const&, double beta, int,
                  double sign) const
      { if (dim > 0) m["Stiffness"] << sign * w2 / (beta * dim); }
    };
    void init_collector(collector& coll) const { coll.init(dim); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m,
                              virtual_lattice_t const& vlat,
                              double beta, double sign,
                              std::vector<int> const& /* spins */,
                              std::vector<OP> const& operators,
                              std::vector<int> const& /* spins_c */,
                              std::vector<FRAGMENT> const& /* fragments */,
                              collector const& coll)
    { coll.commit(m, vlat, beta, operators.size(), sign); }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, virtual_lattice_t const& vlat,
                            double beta, double sign,
                            std::vector<int> const& spins,
                            std::vector<OP> const& operators,
                            std::vector<int>& spins_c)
    {
      if (improved || dim == 0) return;

      std::valarray<double> winding(0., dim);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          if (oi->is_bond()) {
            double s = 1-2*spins_c[vsource(oi->pos(), vlat)];
            alps::coordinate_type const& vr =
              bond_vector_relative[rbond(vlat, oi->pos())];
            for (int i = 0; i < dim; ++i) winding[i] += s * vr[i];
            spins_c[vsource(oi->pos(), vlat)] ^= 1;
            spins_c[vtarget(oi->pos(), vlat)] ^= 1;
          } else {
            spins_c[oi->pos()] ^= 1;
          }
        }
      }

      double w2 = 0;
      for (int i = 0; i < dim; ++i) w2 += power2(winding[i]);
      m["Stiffness"] << sign * w2 / (beta * dim);
    }
  };

  struct evaluator
  {
    static void evaluate(alps::ObservableSet& /* m */,
                         alps::Parameters const& /* params */,
                         alps::ObservableSet const& /* m_in */) {}
  };
};

} // end namespace looper

#endif // LOOPER_STIFFNESS_H
