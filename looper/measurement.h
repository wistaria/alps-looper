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

#ifndef LOOPER_MEASUREMENT_H
#define LOOPER_MEASUREMENT_H

#include "lattice.h"
#include "power.h"
#include "type.h"
#include <alps/alea.h>
#include <boost/call_traits.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <string>
#include <vector>

namespace looper {

//
// helper functions
//

inline
void add_scalar_obs(alps::ObservableSet& m, std::string const& name, bool is_signed = false) {
  if (!m.has(name)) m << make_observable(alps::RealObservable(name), is_signed);
}

inline
void add_vector_obs(alps::ObservableSet& m, std::string const& name, bool is_signed = false) {
  if (!m.has(name)) m << make_observable(alps::RealVectorObservable(name), is_signed);
}

inline
void add_vector_obs(alps::ObservableSet& m, std::string const& name,
  alps::RealVectorObservable::label_type const& label, bool is_signed = false) {
  if (!m.has(name)) m << make_observable(alps::RealVectorObservable(name, label), is_signed);
}


//
// measurement_set
//

namespace { struct null_measurement; }

template<typename M1,
         typename M2 = null_measurement,
         typename M3 = null_measurement,
         typename M4 = null_measurement,
         typename M5 = null_measurement,
         typename M6 = null_measurement,
         typename M7 = null_measurement,
         typename M8 = null_measurement>
struct measurement_set {};


//
// composite_measurement (forward declaration)
//

template<typename M1, typename M2>
struct composite_measurement;


//
// measurement wrapper
//

template<typename ESTIMATOR, typename DUMB_ESTIMATE, typename DUMB_COLLECTOR,
  bool HAS_IMPROVED_ESTIMATOR>
struct improved_estimator_selector_helper {
  typedef DUMB_ESTIMATE estimate;
  typedef DUMB_COLLECTOR collector;
};
template<typename ESTIMATOR, typename DUMB_ESTIMATE, typename DUMB_COLLECTOR>
struct improved_estimator_selector_helper<ESTIMATOR, DUMB_ESTIMATE, DUMB_COLLECTOR, true> {
  typedef typename ESTIMATOR::improved_estimator::estimate estimate;
  typedef typename ESTIMATOR::improved_estimator::collector collector;
};

template<typename MEASUREMENT, typename MC, typename LAT, typename TIME>
struct improved_estimator_selector {
  typedef MEASUREMENT measurement_t;
  typedef typename measurement_t::template estimator<MC, LAT, TIME> estimator_t;

  struct dumb_estimate {
    typedef LAT lattice_t;
    void reset(estimator_t const&) {}
    dumb_estimate& operator+=(dumb_estimate const&) { return *this; }
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
  };

  struct dumb_collector {
    typedef LAT lattice_t;
    void reset(estimator_t const&) {}
    template<typename T>
    dumb_collector& operator+=(T const&) { return *this; }
    template<typename M>
    void commit(M&, lattice_t const&, double /* beta */, double /* sign */,
      double /* nop */) const {}
  };

  typedef improved_estimator_selector_helper<estimator_t, dumb_estimate,
    dumb_collector, boost::is_base_of<has_improved_estimator_tag, measurement_t>::value> helper;
  typedef typename helper::estimate estimate;
  typedef typename helper::collector collector;
};


template<typename ESTIMATOR, typename DUMB_COLLECTOR, bool HAS_NORMAL_ESTIMATOR>
struct normal_estimator_selector_helper {
  typedef DUMB_COLLECTOR collector;
};
template<typename ESTIMATOR, typename DUMB_COLLECTOR>
struct normal_estimator_selector_helper<ESTIMATOR, DUMB_COLLECTOR, true> {
  typedef typename ESTIMATOR::normal_estimator::collector collector;
};

template<typename MEASUREMENT, typename MC, typename LAT, typename TIME>
struct normal_estimator_selector {
  typedef MEASUREMENT measurement_t;
  typedef typename measurement_t::template estimator<MC, LAT, TIME> estimator_t;

  struct dumb_collector {
    typedef LAT lattice_t;
    void reset(estimator_t const&) {}
    dumb_collector& operator+=(dumb_collector const&) { return *this; }
    void begin_s(estimator_t const&, lattice_t const&, double /* t */, int /* s */, int /* c */) {}
    void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
    void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
    void end_s(estimator_t const&, lattice_t const&, double, int, int) {}
    void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) {}
    void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
    void start_bottom(estimator_t const&, lattice_t const&, double, int, int) {}
    void start(estimator_t const&, lattice_t const&, double, int, int) {}
    void stop_top(estimator_t const&, lattice_t const&, double, int, int) {}
    void stop(estimator_t const&, lattice_t const&, double, int, int) {}
    template<typename M>
    void commit(M&, lattice_t const&, double /* beta */, double /* sign */,
      double /* nop */) const {}
  };

  typedef typename normal_estimator_selector_helper<estimator_t, dumb_collector,
    boost::is_base_of<has_normal_estimator_tag, measurement_t>::value>::collector collector;
};

template<typename MEASUREMENT, typename MC, typename LAT, typename TIME,
  bool NEED_IMPROVED_COLLECTOR, bool NEED_NORMALCOLLECTOR>
struct collector_selector {
  typedef typename MEASUREMENT::template estimator<MC, LAT, TIME> estimator_t;
  typedef typename improved_estimator_selector<MEASUREMENT, MC, LAT, TIME>::dumb_collector
    improved_collector;
  typedef typename normal_estimator_selector<MEASUREMENT, MC, LAT, TIME>::dumb_collector
    normal_collector;
  struct collector : public improved_collector, public normal_collector {
    typedef LAT lattice_t;
    collector() : improved_collector(), normal_collector() {}
    void reset(estimator_t const& emt) {
      improved_collector::reset(emt);
      normal_collector::reset(emt);
    }
    collector& operator+=(collector const& rhs) {
      improved_collector::operator+=(rhs);
      normal_collector::operator+=(rhs);
      return *this;
    }
    template<typename T>
    collector& operator+=(T const& rhs) {
      improved_collector::operator+=(rhs);
      return *this;
    }
    template<typename M>
    void commit(M& m, lattice_t const& lat, double beta, double sign, double nop) const {
      improved_collector::commit(m, lat, beta, sign, nop);
      normal_collector::commit(m, lat, beta, sign, nop);
    }
  };
};

template<typename MEASUREMENT, typename MC, typename LAT, typename TIME>
struct collector_selector<MEASUREMENT, MC, LAT, TIME, true, true> {
  typedef typename MEASUREMENT::template estimator<MC, LAT, TIME> estimator_t;
  typedef typename improved_estimator_selector<MEASUREMENT, MC, LAT, TIME>::collector
    improved_collector;
  typedef typename normal_estimator_selector<MEASUREMENT, MC, LAT, TIME>::collector
    normal_collector;
  struct collector : public improved_collector, public normal_collector {
    typedef LAT lattice_t;
    collector() : improved_collector(), normal_collector() {}
    void reset(estimator_t const& emt) {
      improved_collector::reset(emt);
      normal_collector::reset(emt);
    }
    collector& operator+=(collector const& rhs) {
      improved_collector::operator+=(rhs);
      normal_collector::operator+=(rhs);
      return *this;
    }
    template<typename T>
    collector& operator+=(T const& rhs) {
      improved_collector::operator+=(rhs);
      return *this;
    }
    template<typename M>
    void commit(M& m, lattice_t const& lat, double beta, double sign, double nop) const {
      improved_collector::commit(m, lat, beta, sign, nop);
      normal_collector::commit(m, lat, beta, sign, nop);
    }
  };
};

template<typename MEASUREMENT, typename MC, typename LAT, typename TIME>
struct collector_selector<MEASUREMENT, MC, LAT, TIME, true, false> {
  typedef typename MEASUREMENT::template estimator<MC, LAT, TIME> estimator_t;
  typedef typename improved_estimator_selector<MEASUREMENT, MC, LAT, TIME>::collector
    improved_collector;
  typedef typename normal_estimator_selector<MEASUREMENT, MC, LAT, TIME>::dumb_collector
    normal_collector;
  struct collector : public improved_collector, public normal_collector {
    typedef LAT lattice_t;
    collector() : improved_collector(), normal_collector() {}
    void reset(estimator_t const& emt) {
      improved_collector::reset(emt);
      normal_collector::reset(emt);
    }
    collector& operator+=(collector const& rhs) {
      improved_collector::operator+=(rhs);
      normal_collector::operator+=(rhs);
      return *this;
    }
    template<typename T>
    collector& operator+=(T const& rhs) {
      improved_collector::operator+=(rhs);
      return *this;
    }
    template<typename M>
    void commit(M& m, lattice_t const& lat, double beta, double sign, double nop) const {
      improved_collector::commit(m, lat, beta, sign, nop);
      normal_collector::commit(m, lat, beta, sign, nop);
    }
  };
};

template<typename MEASUREMENT, typename MC, typename LAT, typename TIME>
struct collector_selector<MEASUREMENT, MC, LAT, TIME, false, true> {
  typedef typename MEASUREMENT::template estimator<MC, LAT, TIME> estimator_t;
  typedef typename improved_estimator_selector<MEASUREMENT, MC, LAT, TIME>::dumb_collector
    improved_collector;
  typedef typename normal_estimator_selector<MEASUREMENT, MC, LAT, TIME>::collector
    normal_collector;
  struct collector : public improved_collector, public normal_collector {
    typedef LAT lattice_t;
    collector() : improved_collector(), normal_collector() {}
    void reset(estimator_t const& emt) {
      improved_collector::reset(emt);
      normal_collector::reset(emt);
    }
    collector& operator+=(collector const& rhs) {
      improved_collector::operator+=(rhs);
      normal_collector::operator+=(rhs);
      return *this;
    }
    template<typename T>
    collector& operator+=(T const& rhs) {
      improved_collector::operator+=(rhs);
      return *this;
    }
    template<typename M>
    void commit(M& m, lattice_t const& lat, double beta, double sign, double nop) const {
      improved_collector::commit(m, lat, beta, sign, nop);
      normal_collector::commit(m, lat, beta, sign, nop);
    }
  };
};

template<typename MEASUREMENT, typename DUMB_PRE_EVALUATOR, bool HAS_PRE_EVALUATOR>
struct pre_evaluator_selector_helper {
  typedef DUMB_PRE_EVALUATOR pre_evaluator;
};
template<typename MEASUREMENT, typename DUMB_PRE_EVALUATOR>
struct pre_evaluator_selector_helper<MEASUREMENT, DUMB_PRE_EVALUATOR, true> {
  typedef typename MEASUREMENT::pre_evaluator pre_evaluator;
};

template<typename MEASUREMENT>
struct pre_evaluator_selector {
  typedef MEASUREMENT measurement_t;
  struct dumb_pre_evaluator {
    static void pre_evaluate(alps::ObservableSet&, alps::Parameters const&,
      alps::ObservableSet const&) {}
  };
  typedef typename pre_evaluator_selector_helper<measurement_t, dumb_pre_evaluator,
    boost::is_base_of<has_pre_evaluator_tag, MEASUREMENT>::value>::pre_evaluator pre_evaluator;
};

template<typename MEASUREMENT, typename DUMB_EVALUATOR, bool HAS_EVALUATOR>
struct evaluator_selector_helper {
  typedef DUMB_EVALUATOR evaluator;
};
template<typename MEASUREMENT, typename DUMB_EVALUATOR>
struct evaluator_selector_helper<MEASUREMENT, DUMB_EVALUATOR, true> {
  typedef typename MEASUREMENT::evaluator evaluator;
};

template<typename MEASUREMENT>
struct evaluator_selector {
  typedef MEASUREMENT measurement_t;
  struct dumb_evaluator {
    static void evaluate(alps::ObservableSet&, alps::Parameters const&,
      alps::ObservableSet const&) {}
  };
  typedef typename evaluator_selector_helper<measurement_t, dumb_evaluator,
    boost::is_base_of<has_evaluator_tag, MEASUREMENT>::value>::evaluator evaluator;
};

template<typename MEASUREMENT>
struct measurement_wrapper {
  typedef MEASUREMENT measurement_t;
  template<typename MC, typename LAT, typename TIME>
  struct estimator : public measurement_t::template estimator<MC, LAT, TIME> {
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    struct improved_estimator {
      typedef typename improved_estimator_selector<measurement_t, mc_type, lattice_t, time_t>::
        estimate estimate;
      typedef typename collector_selector<measurement_t, mc_type, lattice_t, time_t,
        boost::is_base_of<has_improved_estimator_tag, MEASUREMENT>::value,
        (!boost::is_base_of<has_improved_estimator_tag, MEASUREMENT>::value &&
         boost::is_base_of<has_normal_estimator_tag, MEASUREMENT>::value)>::collector collector;
    };
    struct normal_estimator {
      typedef typename collector_selector<measurement_t, mc_type, lattice_t, time_t,
        false, boost::is_base_of<has_normal_estimator_tag, MEASUREMENT>::value>::collector
        collector;
    };
    // inputs:
    //   has_improved_estimator        1    1    1    1    0    0    0    0
    //   has_normal_estimator          1    1    0    0    1    1    0    0
    //   enable_improved_estimator     1    0    1    0    1    0    1    0
    // results:
    //   improved collector            I    d    I    d    d    d    d    d
    //   normal collector              d    N    d    d    N    N    d    d
  };
  typedef typename pre_evaluator_selector<measurement_t>::pre_evaluator pre_evaluator;
  typedef typename evaluator_selector<measurement_t>::evaluator evaluator;
};


//
// measurement traits
//

template<typename MEASUREMENT_SET>
struct measurement;

template<typename M1>
struct measurement<measurement_set<M1, null_measurement, null_measurement, null_measurement,
  null_measurement, null_measurement, null_measurement, null_measurement> > {
  typedef measurement_wrapper<M1> type;
};

template<typename M1, typename M2>
struct measurement<measurement_set<M1, M2, null_measurement, null_measurement, null_measurement,
  null_measurement, null_measurement, null_measurement> > {
  typedef composite_measurement<measurement_wrapper<M1>, measurement_wrapper<M2> > type;
};

template<typename M1, typename M2, typename M3>
struct measurement<measurement_set<M1, M2, M3, null_measurement, null_measurement,
  null_measurement, null_measurement, null_measurement> > {
  typedef composite_measurement<composite_measurement<measurement_wrapper<M1>,
    measurement_wrapper<M2> >, measurement_wrapper<M3> > type;
};

template<typename M1, typename M2, typename M3, typename M4>
struct measurement<measurement_set<M1, M2, M3, M4, null_measurement, null_measurement,
  null_measurement, null_measurement> > {
  typedef composite_measurement<composite_measurement<composite_measurement<
    measurement_wrapper<M1>, measurement_wrapper<M2> >, measurement_wrapper<M3> >,
    measurement_wrapper<M4> > type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5>
struct measurement<measurement_set<M1, M2, M3, M4, M5, null_measurement, null_measurement,
  null_measurement> > {
  typedef composite_measurement<composite_measurement<composite_measurement<composite_measurement<
    measurement_wrapper<M1>, measurement_wrapper<M2> >, measurement_wrapper<M3> >,
    measurement_wrapper<M4> >, measurement_wrapper<M5> > type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5, typename M6>
struct measurement<measurement_set<M1, M2, M3, M4, M5, M6, null_measurement, null_measurement> > {
  typedef composite_measurement<composite_measurement<composite_measurement<composite_measurement<
    composite_measurement<measurement_wrapper<M1>, measurement_wrapper<M2> >,
    measurement_wrapper<M3> >, measurement_wrapper<M4> >, measurement_wrapper<M5> >,
    measurement_wrapper<M6> > type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7>
struct measurement<measurement_set<M1, M2, M3, M4, M5, M6, M7, null_measurement> > {
  typedef composite_measurement<composite_measurement<composite_measurement<composite_measurement<
    composite_measurement<composite_measurement<measurement_wrapper<M1>, measurement_wrapper<M2> >,
    measurement_wrapper<M3> >, measurement_wrapper<M4> >, measurement_wrapper<M5> >,
    measurement_wrapper<M6> >, measurement_wrapper<M7> > type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7,
  typename M8>
struct measurement<measurement_set<M1, M2, M3, M4, M5, M6, M7, M8> > {
  typedef composite_measurement<composite_measurement<composite_measurement<composite_measurement<
    composite_measurement<composite_measurement<composite_measurement<measurement_wrapper<M1>,
    measurement_wrapper<M2> >, measurement_wrapper<M3> >, measurement_wrapper<M4> >,
    measurement_wrapper<M5> >, measurement_wrapper<M6> >, measurement_wrapper<M7> >,
    measurement_wrapper<M8> > type;
};


//
// selectors
//

template<typename MEASUREMENT_SET, typename MC, typename LAT, typename TIME>
struct estimator {
  typedef typename measurement<MEASUREMENT_SET>::type measurement_t;
  typedef typename measurement_t::template estimator<MC, LAT, TIME> type;
};


//
// struct base_measurement
//

  struct base_measurement : public has_improved_estimator_tag, public has_normal_estimator_tag {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef estimator<MC, LAT, TIME> estimator_t;
    typedef LAT lattice_t;
    void initialize(alps::Parameters const&, lattice_t const&, bool, bool) {}
    template<typename M>
    void init_observables(M&, bool, bool) {}

    struct collector_common {
      double nop_; // total number of operators
      double nc_; // total number of (closed) clusters
      std::pair<int, int> range_; // configuration range (for parallel QMC)
      unsigned int noc_; // number of open clusters (for parallel QMC)
      collector_common() : nop_(0), nc_(0), range_(std::make_pair(1, 0)), noc_(0) {}
      void reset(estimator_t const&) {
        nop_ = 0; nc_= 0; range_ = std::make_pair(1, 0);  noc_ = 0;
      }
      collector_common& operator+=(collector_common const& rhs) {
        nop_ += rhs.nop_;
        nc_ += rhs.nc_;
        range_ = std::make_pair(std::min(range_.first, rhs.range_.first),
                                std::max(range_.second, rhs.range_.second));
        return *this;
      }
      void set_num_clusters(unsigned int n) { nc_ = n; }
      void inc_num_clusters(unsigned int n) { nc_ += n; }
      double num_clusters() const { return nc_; }
      void set_num_operators(unsigned int n) { nop_ = n; }
      double num_operators() const { return nop_; }
      // for parallel QMC
      void set_num_open_clusters(unsigned int n) { noc_ = n; }
      unsigned int num_open_clusters() const { return noc_; }
      void clear_range() { range_ = std::make_pair(1, 0); }
      void set_range(int pos) { range_ = std::make_pair(pos, pos); }
      std::pair<int, int> const& range() const { return range_; }
      bool empty() const { return range_.first > range_.second; }
      template<typename M>
      void commit(M&, lattice_t const&, double, double, double) const {}
    };

    struct improved_estimator {
      struct estimate {
        bool to_flip;
        void reset(estimator_t const&) {}
        estimate& operator+=(estimate const& rhs) {
          to_flip ^= rhs.to_flip;
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
      };
      struct collector : public collector_common {
        typedef collector_common base_type;
        collector() : base_type() {}
        collector& operator+=(estimate const&) { return *this; }
        collector& operator+=(collector const& rhs) {
          base_type::operator+=(rhs);
          return *this;
        }
      };
    };

    struct normal_estimator {
      struct collector : public collector_common {
        typedef collector_common base_type;
        collector() : base_type() {}
        collector& operator+=(collector const& rhs) {
          base_type::operator+=(rhs);
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
      };
    };
  };
};


//
// accumulator
//

template<typename ESTIMATOR, typename ENABLE_IMPROVED_ESTIMATOR>
struct normal_accumulator {
  typedef ESTIMATOR estimator_t;
  typedef typename estimator_t::normal_estimator::collector collector_t;
  typedef typename estimator_t::lattice_t lattice_t;
  typedef typename boost::call_traits<typename estimator_t::time_t>::param_type time_pt;
  normal_accumulator(collector_t& cl, lattice_t const& lt, estimator_t const& emt) :
    coll(cl), lat(lt), estimator(emt) {}
  void begin_s(time_pt t, int s, int c) {
    coll.begin_s(estimator, lat, t, s, c);
  }
  void begin_b(time_pt t, int b, int s0, int s1, int c0, int c1) {
    coll.begin_bs(estimator, lat, t, b, s0, c0);
    coll.begin_bt(estimator, lat, t, b, s1, c1);
  }
  void end_s(time_pt t, int s, int c) {
    coll.end_s(estimator, lat, t, s, c);
  }
  void end_b(time_pt t, int b, int s0, int s1, int c0, int c1) {
    coll.end_bs(estimator, lat, t, b, s0, c0);
    coll.end_bt(estimator, lat, t, b, s1, c1);
  }
  void start_bottom(time_pt t, int s, int c) {
    coll.start_bottom(estimator, lat, t, s, c);
  }
  void start(time_pt t, int s, int c) {
    coll.start(estimator, lat, t, s, c);
  }
  void stop(time_pt t, int s, int c) {
    coll.stop(estimator, lat, t, s, c);
  }
  void stop_top(time_pt t, int s, int c) {
    coll.stop_top(estimator, lat, t, s, c);
  }
  collector_t& coll;
  lattice_t const& lat;
  estimator_t const& estimator;
};

template<typename ESTIMATOR>
struct normal_accumulator<ESTIMATOR, boost::mpl::true_> {
  typedef ESTIMATOR estimator_t;
  typedef typename estimator_t::improved_estimator::collector collector_t;
  typedef typename estimator_t::lattice_t lattice_t;
  typedef typename boost::call_traits<typename estimator_t::time_t>::param_type time_pt;
  normal_accumulator(collector_t&, lattice_t const&, estimator_t const&) {}
  void begin_s(time_pt, int, int) {}
  void begin_b(time_pt, int, int, int, int, int) {}
  void end_s(time_pt, int, int) {}
  void end_b(time_pt, int, int, int, int, int) {}
  void start_bottom(time_pt, int, int) {}
  void start(time_pt, int, int) {}
  void stop(time_pt, int, int) {}
  void stop_top(time_pt, int, int) {}
};


template<typename ESTIMATOR, typename FRAGMENT, typename ESTIMATE,
  typename ENABLE_IMPROVED_ESTIMATOR>
struct improved_accumulator {
  typedef ESTIMATOR estimator_t;
  typedef FRAGMENT fragment_t;
  typedef ESTIMATE estimate_t;
  typedef typename estimator_t::lattice_t lattice_t;
  typedef typename boost::call_traits<typename estimator_t::time_t>::param_type time_pt;
  improved_accumulator(std::vector<estimate_t>&, lattice_t const&, estimator_t const&,
    std::vector<fragment_t> const&) {}
  void begin_s(int, time_pt, int, int) {}
  void begin_b(int, int, time_pt, int, int, int, int, int) {}
  void end_s(int, time_pt, int, int) {}
  void end_b(int, int, time_pt, int, int, int, int, int) {}
  void start_bottom(int, time_pt, int, int) {}
  void start(int, time_pt, int, int) {}
  void stop(int, time_pt, int, int) {}
  void stop_top(int, time_pt, int, int) {}
};

template<typename ESTIMATOR, typename FRAGMENT, typename ESTIMATE>
struct improved_accumulator<ESTIMATOR, FRAGMENT, ESTIMATE, boost::mpl::true_> {
  typedef ESTIMATOR estimator_t;
  typedef FRAGMENT fragment_t;
  typedef ESTIMATE estimate_t;
  typedef typename estimator_t::lattice_t lattice_t;
  typedef typename boost::call_traits<typename estimator_t::time_t>::param_type time_pt;
  improved_accumulator(std::vector<estimate_t>& es, lattice_t const& lt,
    estimator_t const& emt, std::vector<fragment_t> const& fr) :
    estimates(es), lat(lt), estimator(emt), fragments(fr) {}
  void begin_s(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].begin_s(estimator, lat, t, s, c);
  }
  void begin_b(int p0, int p1, time_pt t, int b, int s0, int s1, int c0, int c1) {
    estimates[fragments[p0].id()].begin_bs(estimator, lat, t, b, s0, c0);
    estimates[fragments[p1].id()].begin_bt(estimator, lat, t, b, s1, c1);
  }
  void end_s(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].end_s(estimator, lat, t, s, c);
  }
  void end_b(int p0, int p1, time_pt t, int b, int s0, int s1, int c0, int c1) {
    estimates[fragments[p0].id()].end_bs(estimator, lat, t, b, s0, c0);
    estimates[fragments[p1].id()].end_bt(estimator, lat, t, b, s1, c1);
  }
  void start_bottom(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].start_bottom(estimator, lat, t, s, c);
  }
  void start(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].start(estimator, lat, t, s, c);
  }
  void stop(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].stop(estimator, lat, t, s, c);
  }
  void stop_top(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].stop_top(estimator, lat, t, s, c);
  }
  std::vector<estimate_t>& estimates;
  lattice_t const& lat;
  estimator_t const& estimator;
  std::vector<fragment_t> const& fragments;
};


//
// composite_measurement
//

template<typename MEASUREMENT1, typename MEASUREMENT2>
struct composite_measurement {
  typedef MEASUREMENT1 measurement1;
  typedef MEASUREMENT2 measurement2;

  template<typename MC, typename LAT, typename TIME>
  struct estimator : public measurement1::template estimator<MC, LAT, TIME>,
    public measurement2::template estimator<MC, LAT, TIME> {
    typedef MC mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;

    typedef typename measurement1::template estimator<MC, LAT, TIME> estimator1;
    typedef typename measurement2::template estimator<MC, LAT, TIME> estimator2;

    typedef typename boost::call_traits<time_t>::param_type time_pt;

    void initialize(alps::Parameters const& params, lattice_t const& lat, bool is_signed,
      bool enable_improved_estimator) {
      estimator1::initialize(params, lat, is_signed, enable_improved_estimator);
      estimator2::initialize(params, lat, is_signed, enable_improved_estimator);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed, bool enable_improved_estimator) {
      estimator1::init_observables(m, is_signed, enable_improved_estimator);
      estimator2::init_observables(m, is_signed, enable_improved_estimator);
    }

    struct improved_estimator {
      typedef typename base_measurement::estimator<MC, LAT, TIME>::improved_estimator::estimate
        minimal_estimate;
      typedef typename estimator1::improved_estimator::estimate estimate1;
      typedef typename estimator2::improved_estimator::estimate estimate2;
      struct estimate : public estimate1, public estimate2 {
        estimate() : estimate1(), estimate2() {}
        void reset(estimator_t const& emt) {
          estimate1::reset(emt);
          estimate2::reset(emt);
        }
        estimate& operator+=(estimate const& rhs) {
          estimate1::operator+=(rhs);
          estimate2::operator+=(rhs);
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          estimate1::begin_s(emt, lat, t, s, c);
          estimate2::begin_s(emt, lat, t, s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s,
          int c) {
          estimate1::begin_bs(emt, lat, t, b, s, c);
          estimate2::begin_bs(emt, lat, t, b, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s,
          int c) {
          estimate1::begin_bt(emt, lat, t, b, s, c);
          estimate2::begin_bt(emt, lat, t, b, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          estimate1::end_s(emt, lat, t, s, c);
          estimate2::end_s(emt, lat, t, s, c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
          estimate1::end_bs(emt, lat, t, b, s, c);
          estimate2::end_bs(emt, lat, t, b, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
          estimate1::end_bt(emt, lat, t, b, s, c);
          estimate2::end_bt(emt, lat, t, b, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          estimate1::start_bottom(emt, lat, t, s, c);
          estimate2::start_bottom(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          estimate1::start(emt, lat, t, s, c);
          estimate2::start(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          estimate1::stop(emt, lat, t, s, c);
          estimate2::stop(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          estimate1::stop_top(emt, lat, t, s, c);
          estimate2::stop_top(emt, lat, t, s, c);
        }
      };

      typedef typename estimator1::improved_estimator::collector collector1;
      typedef typename estimator2::improved_estimator::collector collector2;
      struct collector : public collector1, public collector2 {
        collector() : collector1(), collector2() {}
        void reset(estimator_t const& emt) {
          collector1::reset(emt);
          collector2::reset(emt);
        }
        collector& operator+=(collector const& rhs) {
          collector1::operator+=(rhs);
          collector2::operator+=(rhs);
          return *this;
        }
        collector& operator+=(estimate const& rhs) {
          collector1::operator+=(rhs);
          collector2::operator+=(rhs);
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::begin_s(emt, lat, t, s, c);
          collector2::begin_s(emt, lat, t, s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s,
          int c) {
          collector1::begin_bs(emt, lat, t, b, s, c);
          collector2::begin_bs(emt, lat, t, b, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s,
          int c) {
          collector1::begin_bt(emt, lat, t, b, s, c);
          collector2::begin_bt(emt, lat, t, b, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::end_s(emt, lat, t, s, c);
          collector2::end_s(emt, lat, t, s, c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
          collector1::end_bs(emt, lat, t, b, s, c);
          collector2::end_bs(emt, lat, t, b, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
          collector1::end_bt(emt, lat, t, b, s, c);
          collector2::end_bt(emt, lat, t, b, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::start_bottom(emt, lat, t, s, c);
          collector2::start_bottom(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::start(emt, lat, t, s, c);
          collector2::start(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::stop(emt, lat, t, s, c);
          collector2::stop(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::stop_top(emt, lat, t, s, c);
          collector2::stop_top(emt, lat, t, s, c);
        }
        template<typename M>
        void commit(M& m, lattice_t const& lat, double beta, double sign, double nop) const {
          collector1::commit(m, lat, beta, sign, nop);
          collector2::commit(m, lat, beta, sign, nop);
        }
      };
    };

    struct normal_estimator {
      typedef typename estimator1::normal_estimator::collector collector1;
      typedef typename estimator2::normal_estimator::collector collector2;
      struct collector : public collector1, public collector2 {
        collector() : collector1(), collector2() {}
        void reset(estimator_t const& emt) {
          collector1::reset(emt);
          collector2::reset(emt);
        }
        collector& operator+=(collector const& rhs) {
          collector1::operator+=(rhs);
          collector2::operator+=(rhs);
          return *this;
        }
        template<typename T>
        collector& operator+=(T const& rhs) {
          collector1::operator+=(rhs);
          collector2::operator+=(rhs);
          return *this;
        }
        void begin_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::begin_s(emt, lat, t, s, c);
          collector2::begin_s(emt, lat, t, s, c);
        }
        void begin_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s,
          int c) {
          collector1::begin_bs(emt, lat, t, b, s, c);
          collector2::begin_bs(emt, lat, t, b, s, c);
        }
        void begin_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s,
          int c) {
          collector1::begin_bt(emt, lat, t, b, s, c);
          collector2::begin_bt(emt, lat, t, b, s, c);
        }
        void end_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::end_s(emt, lat, t, s, c);
          collector2::end_s(emt, lat, t, s, c);
        }
        void end_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
          collector1::end_bs(emt, lat, t, b, s, c);
          collector2::end_bs(emt, lat, t, b, s, c);
        }
        void end_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
          collector1::end_bt(emt, lat, t, b, s, c);
          collector2::end_bt(emt, lat, t, b, s, c);
        }
        void start_bottom(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::start_bottom(emt, lat, t, s, c);
          collector2::start_bottom(emt, lat, t, s, c);
        }
        void start(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::start(emt, lat, t, s, c);
          collector2::start(emt, lat, t, s, c);
        }
        void stop(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::stop(emt, lat, t, s, c);
          collector2::stop(emt, lat, t, s, c);
        }
        void stop_top(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
          collector1::stop_top(emt, lat, t, s, c);
          collector2::stop_top(emt, lat, t, s, c);
        }
        template<typename M>
        void commit(M& m, lattice_t const& lat, double beta, double sign, double nop) const {
          collector1::commit(m, lat, beta, sign, nop);
          collector2::commit(m, lat, beta, sign, nop);
        }
      };
    };
  };

  typedef typename measurement1::pre_evaluator pre_evaluator1;
  typedef typename measurement2::pre_evaluator pre_evaluator2;
  struct pre_evaluator : public pre_evaluator1, public pre_evaluator2 {
    static void pre_evaluate(alps::ObservableSet& m, alps::Parameters const& params,
      alps::ObservableSet const& m_in) {
      pre_evaluator1::pre_evaluate(m, params, m_in);
      pre_evaluator2::pre_evaluate(m, params, m_in);
    }
  };

  typedef typename measurement1::evaluator evaluator1;
  typedef typename measurement2::evaluator evaluator2;
  struct evaluator : public evaluator1, public evaluator2 {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& params,
      alps::ObservableSet const& m_in) {
      evaluator1::evaluate(m, params, m_in);
      evaluator2::evaluate(m, params, m_in);
    }
  };
};

} // end namespace looper

#endif // LOOPER_MEASUREMENT_H
