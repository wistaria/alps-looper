/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

// $Id: percolation_impl.h 660 2004-03-04 10:58:14Z wistaria $
// percolation_impl.h - implementation of worker for percolation simulation

#ifndef PERCOLATION_IMPL_H
#define PERCOLATION_IMPL_H

#include <looper/copyright.h>
#include <looper/union_find.h>

#include <alps/alea.h>
#include <alps/scheduler.h>
#include <boost/graph/graph_traits.hpp>
#include <utility>

struct percolation
{
  struct node_base {
    bool occupied;
    void reset() { occupied = true; }
    node_base& operator+=(const node_base&) { return *this; }
  };
  typedef looper::union_find::node<node_base> node_type;
  typedef std::vector<node_type> vector_type;


  template<class RNG, bool BOND_P> struct initializer;
  
  template<class RNG>
  struct initializer<RNG, false>
  {
    initializer(double p, RNG& rng) : p_(p), rng_(rng) {}
    template<class T>
    void operator()(T& a) const { a.reset(); a.occupied = (rng_() < p_); }
    double p_;
    mutable RNG& rng_;
  };

  template<class RNG>
  struct initializer<RNG, true>
  {
    template<class T>
    void operator()(T& a) const { a.reset(); }
  };


  template<class G, class V, class RNG, bool BOND_P> struct unifier;
  
  template<class G, class V, class RNG>
  struct unifier<G, V, RNG, false>
  {
    unifier(const G& g, V& v) : g_(g), v_(v) {}
    template<class T>
    void operator()(const T& a) const {
      typename V::iterator s_itr = v_.begin() + boost::source(a, g_);
      typename V::iterator t_itr = v_.begin() + boost::target(a, g_);
      if (s_itr->occupied && t_itr->occupied)
	looper::union_find::unify(*s_itr, *t_itr);
    }
    const G& g_;
    mutable V& v_;
  };

  template<class G, class V, class RNG>
  struct unifier<G, V, RNG, true>
  {
    unifier(const G& g, V& v, double p, RNG& rng) :
      g_(g), v_(v), p_(p), rng_(rng) {}
    template<class T>
    void operator()(const T& a) const {
      if (rng_() < p_) {
	typename V::iterator s_itr = v_.begin() + boost::source(a, g_);
	typename V::iterator t_itr = v_.begin() + boost::target(a, g_);
	looper::union_find::unify(*s_itr, *t_itr);
      }
    }
    const G& g_;
    mutable V& v_;
    double p_;
    mutable RNG& rng_;
  };


  class worker_base
  {
  public:
    typedef alps::SimpleRealObservable measurement_type;
    
    template<class G>
    worker_base(const G& g, double p, bool bp) :
      bond_p_(bp), vertices_(boost::num_vertices(g)), prob_(p) {}
    
    template<class G>
    worker_base(const G& g, double p, bool bp, alps::ObservableSet& m) :
      bond_p_(bp), vertices_(boost::num_vertices(g)), prob_(p)
    {
      // setup measurements
      m << measurement_type("Number of Clusters");
      m << measurement_type("Percolation Probability");
      m << measurement_type("Disconnected Susceptiblity");
      m << measurement_type("Connected Susceptiblity");
      m.reset(true);
    }
    
    template<class G, class RNG>
    void step(const G& g, RNG& random_01, alps::ObservableSet& m)
    {
      // setup sample
      if (!bond_p_)
	std::for_each(vertices_.begin(), vertices_.end(),
		      initializer<RNG, false>(prob_, random_01));
      else
	std::for_each(vertices_.begin(), vertices_.end(),
		      initializer<RNG, true>());
      
      // union_find
      if (!bond_p_)
 	std::for_each(boost::edges(g).first, boost::edges(g).second,
 		      unifier<G, vector_type, RNG, false>(g, vertices_));
      else
 	std::for_each(boost::edges(g).first, boost::edges(g).second,
 		      unifier<G, vector_type, RNG, true>(g, vertices_, prob_,
 							 random_01));
      
      // measurement
      int num_clusters = 0;
      double max_weight = 0.;
      double sc = 0.;
      vector_type::iterator itr_end = vertices_.end();
      for (vector_type::iterator itr = vertices_.begin();
	   itr != itr_end; ++itr) {
	if (itr->occupied && itr->is_root()) {
	  ++num_clusters;
	  double w = itr->weight();
	  max_weight = std::max(max_weight, w);
	  sc += w * w;
	}
      }
      
      double nv = double(boost::num_vertices(g));
      m.template get<measurement_type>("Number of Clusters") <<
	double(num_clusters) / nv;
      m.template get<measurement_type>("Percolation Probability") <<
	max_weight / nv;
      m.template get<measurement_type>("Disconnected Susceptiblity") <<
	sc / nv;
      m.template get<measurement_type>("Connected Susceptiblity") <<
	(sc - max_weight * max_weight) / nv;
    }

    void output_results(std::ostream& os, alps::ObservableSet& m) const {
      os << m.get<measurement_type>("Number of Clusters").mean() << ' '
	 << m.get<measurement_type>("Number of Clusters").error() << ' '
	 << m.get<measurement_type>("Percolation Probability").mean() << ' '
	 << m.get<measurement_type>("Percolation Probability").error() << ' '
	 << m.get<measurement_type>("Disconnected Susceptiblity").mean() << ' '
	 << m.get<measurement_type>("Disconnected Susceptiblity").error() << ' '
	 << m.get<measurement_type>("Connected Susceptiblity").mean() << ' '
	 << m.get<measurement_type>("Connected Susceptiblity").error();
    }
  private:
    bool bond_p_;
    vector_type vertices_;
    double prob_;
  };


  class worker : public alps::scheduler::LatticeMCRun<>
  {
  public:
    worker(const alps::ProcessList& w, const alps::Parameters& p, int n) :
      alps::scheduler::LatticeMCRun<>(w, p, n),
      prob_(double(p["P"])),
      samples_(int(p["SAMPLES"])), samples_done_(0),
      wb_(graph(), prob_, (p["TYPE"] == "bond"), measurements) {}
    virtual ~worker() {}
    
    virtual void dostep() {
      ++samples_done_;
      wb_.step(graph(), random_01, measurements);
    }
    bool is_thermalized() const { return true; }
    virtual double work_done() const {
      return double(samples_done_) / samples_;
    }
    
    virtual void save(alps::ODump& od) const {
      od << prob_ << samples_ << samples_done_;
    }
    virtual void load(alps::IDump& id) {
      id >> prob_ >> samples_ >> samples_done_;
      if (where.empty()) measurements.compact();
    }

  private:
    double prob_;
    unsigned int samples_;
    unsigned int samples_done_;
    worker_base wb_;
  };
  
  class factory : public alps::scheduler::Factory
  {
    alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
      const boost::filesystem::path& fn) const {
      return new alps::scheduler::MCSimulation(w, fn);
    }
    alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
      const boost::filesystem::path& fn, const alps::Parameters&) const {
      return new alps::scheduler::MCSimulation(w, fn);
    }

    worker* make_worker(const alps::ProcessList& w,
			const alps::Parameters& p, int n) const {
      return new worker(w, p, n);
    }

    void print_copyright(std::ostream& os) const {
      looper::print_copyright(os);
    }
  };
};

#endif // PERCOLATION_IMPL_H
