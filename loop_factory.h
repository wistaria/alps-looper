/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/evaluator.h>
#ifdef HAVE_PARAPACK
# include <parapack/exchange.h>
# include <parapack/worker.h>
#endif
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/throw_exception.hpp>
#include <map>
#include <stdexcept>

template<typename WORKER>
class alps_worker_wrapper : public alps::scheduler::MCRun {
public:
  alps_worker_wrapper(alps::ProcessList const& w, alps::Parameters const& p, int n) :
    alps::scheduler::MCRun(w, p, n),  worker_(p, measurements) {}
  virtual ~alps_worker_wrapper() {}
  virtual void dostep() { worker_.run(*engine_ptr, measurements); }
  virtual void load(alps::IDump& dump) { worker_.load(dump); }
  virtual void save(alps::ODump& dump) const { worker_.save(dump); }

  virtual bool is_thermalized() const { return worker_.is_thermalized(); }
  virtual double work_done() const { return worker_.progress(); }
private:
  WORKER worker_;
};

#ifdef HAVE_PARAPACK

template<typename WORKER>
class parapack_worker_wrapper : public alps::parapack::mc_worker {
public:
  parapack_worker_wrapper(alps::Parameters const& p, alps::ObservableSet& obs) :
    alps::parapack::mc_worker(p),  worker_(p, obs) {}
  virtual ~parapack_worker_wrapper() {}
  virtual bool is_thermalized() const { return worker_.is_thermalized(); }
  virtual double progress() const { return worker_.progress(); }
  virtual void run(alps::ObservableSet& obs) { worker_.run(*engine_ptr, obs); }
  void set_beta(double beta) { worker_.set_beta(beta); }
  double g_weight() const { return worker_.g_weight(); }
  double lambda(double beta) const { return worker_.lambda(beta); }
  virtual void load(alps::IDump& dump) { worker_.load(dump); }
  virtual void save(alps::ODump& dump) const { worker_.save(dump); }
private:
  WORKER worker_;
};

#endif // HAVE_PARAPACK

class abstract_worker_creator {
public:
  virtual ~abstract_worker_creator() {}
  virtual alps::scheduler::MCRun* create(alps::ProcessList const& w,
    alps::Parameters const& p, int n) const = 0;
#ifdef HAVE_PARAPACK
  virtual alps::parapack::abstract_worker* create(const alps::Parameters& p,
    std::vector<alps::ObservableSet>& obs) const = 0;
#endif // HAVE_PARAPACK
};

template <typename WORKER>
class worker_creator : public abstract_worker_creator {
public:
  typedef WORKER worker_type;
  virtual ~worker_creator() {}
  alps::scheduler::MCRun* create(alps::ProcessList const& w, alps::Parameters const& p,
    int n) const {
    return new alps_worker_wrapper<worker_type>(w, p, n);
  }
#ifdef HAVE_PARAPACK
  virtual alps::parapack::abstract_worker* create(const alps::Parameters& p,
    std::vector<alps::ObservableSet>& obs) const {
    if (p.defined("EXCHANGE_MONTE_CARLO") && static_cast<bool>(p["EXCHANGE_MONTE_CARLO"])) {
      return new alps::parapack::single_exchange_worker<parapack_worker_wrapper<worker_type> >
        (p, obs);
    } else {
      obs.resize(1);
      return new parapack_worker_wrapper<worker_type>(p, obs[0]);
    }
  }
#endif // HAVE_PARAPACK
};

class abstract_evaluator_creator {
public:
  virtual ~abstract_evaluator_creator() {}
  virtual looper::abstract_evaluator* create() const = 0;
};

template <typename EVALUATOR>
class evaluator_creator : public abstract_evaluator_creator {
public:
  typedef EVALUATOR evaluator_type;
  virtual ~evaluator_creator() {}
  looper::abstract_evaluator* create() const { return new evaluator_type(); }
};


class loop_factory :
  public alps::scheduler::Factory,
#ifdef HAVE_PARAPACK
  public alps::parapack::abstract_single_worker_factory,
#endif
  private boost::noncopyable {
public:
  alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
    const boost::filesystem::path& fn) const;
  alps::scheduler::MCSimulation* make_task(const alps::ProcessList& w,
    const boost::filesystem::path& fn, const alps::Parameters&) const;
  alps::scheduler::MCSimulation* make_task(const alps::ProcessList&,
    alps::Parameters const&) const {
    return 0;
  }
  alps::scheduler::MCRun* make_worker(alps::ProcessList const& w, alps::Parameters const& p,
    int n) const;
#ifdef HAVE_PARAPACK
  alps::parapack::abstract_worker* make_worker(alps::Parameters const& params,
    std::vector<alps::ObservableSet>& obs) const;
#endif // HAVE_PARAPACK

  looper::abstract_evaluator* make_evaluator(alps::Parameters const& p) const;
  void pre_evaluate(std::vector<alps::ObservableSet>& obs, alps::Parameters const& p) const;
  void evaluate(std::vector<alps::ObservableSet>& obs, alps::Parameters const& p) const;

  void print_copyright(std::ostream& os) const;

  template<typename WORKER>
  bool register_worker(std::string const& name) {
    bool isnew = (worker_creators_.find(name) == worker_creators_.end());
    worker_creators_[name] = worker_pointer_type(new worker_creator<WORKER>());
    return isnew;
  }

  bool unregister_worker(std::string const& name) {
    worker_map_type::iterator itr = worker_creators_.find(name);
    if (itr == worker_creators_.end()) return false;
    worker_creators_.erase(itr);
    return true;
  }

  template<typename EVALUATOR>
  bool register_evaluator(std::string const& name) {
    bool isnew = (evaluator_creators_.find(name) == evaluator_creators_.end());
    evaluator_creators_[name] = evaluator_pointer_type(new evaluator_creator<EVALUATOR>());
    return isnew;
  }

  bool unregister_evaluator(std::string const& name) {
    evaluator_map_type::iterator itr = evaluator_creators_.find(name);
    if (itr == evaluator_creators_.end()) return false;
    evaluator_creators_.erase(itr);
    return true;
  }

  static loop_factory* instance() {
    if (!instance_) instance_ = new loop_factory;
    return instance_;
  }

private:
  loop_factory() {}
  ~loop_factory() {}

  static loop_factory* instance_;

  typedef boost::shared_ptr<abstract_worker_creator> worker_pointer_type;
  typedef std::map<std::string, worker_pointer_type> worker_map_type;
  worker_map_type worker_creators_;

  typedef boost::shared_ptr<abstract_evaluator_creator> evaluator_pointer_type;
  typedef std::map<std::string, evaluator_pointer_type> evaluator_map_type;
  evaluator_map_type evaluator_creators_;
};
