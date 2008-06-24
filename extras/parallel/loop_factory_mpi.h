/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

#include <parapack/parallel.h>
// #include <parapack/exchange_mpi.h>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/throw_exception.hpp>
#include <map>
#include <stdexcept>

template<typename WORKER>
class parallel_worker : public alps::parapack::mc_worker {
public:
  typedef double weight_parameter_type;
  parallel_worker(alps::communicator_helper const& comm, alps::Parameters const& p,
    alps::ObservableSet& obs) : alps::parapack::mc_worker(p),  worker_(comm, p, obs) {}
  virtual ~parallel_worker() {}
  virtual bool is_thermalized() const { return worker_.is_thermalized(); }
  virtual double progress() const { return worker_.progress(); }
  virtual void run(alps::ObservableSet& obs) { worker_.run(*engine_ptr, obs); }
  void set_beta(double beta) { worker_.set_beta(beta); }
  double weight_parameter() const { return worker_.weight_parameter(); }
  static double log_weight(double gw, double beta) { return WORKER::log_weight(gw, beta); }
  virtual void load(alps::IDump& dump) { worker_.load(dump); }
  virtual void save(alps::ODump& dump) const { worker_.save(dump); }
private:
  WORKER worker_;
};

class abstract_parallel_worker_creator {
public:
  virtual ~abstract_parallel_worker_creator() {}
  virtual alps::parapack::abstract_worker* create(alps::communicator_helper const& comm,
    const alps::Parameters& p, std::vector<alps::ObservableSet>& obs) const = 0;
};

template <typename WORKER>
class parallel_worker_creator : public abstract_parallel_worker_creator {
public:
  typedef WORKER worker_type;
  virtual ~parallel_worker_creator() {}
  virtual alps::parapack::abstract_worker* create(alps::communicator_helper const& comm,
    const alps::Parameters& p, std::vector<alps::ObservableSet>& obs) const {
//     if (p.defined("EXCHANGE_MONTE_CARLO") && static_cast<bool>(p["EXCHANGE_MONTE_CARLO"])) {
//       return new alps::parapack::parallel_exchange_worker<parallel_worker<worker_type> >
//         (comm, p, obs);
//     } else {
      obs.resize(1);
      return new parallel_worker<worker_type>(comm, p, obs[0]);
//     }
  }
};


class parallel_worker_factory :
  public alps::parapack::abstract_parallel_worker_factory, private boost::noncopyable {
public:
  alps::parapack::abstract_worker* make_worker(alps::communicator_helper const& comm,
    alps::Parameters const& params, std::vector<alps::ObservableSet>& obs) const;

  std::string version() const;
  void print_copyright(std::ostream& os) const;

  void pre_evaluate(std::vector<alps::ObservableSet>&, alps::Parameters const&) const {}
  void evaluate(std::vector<alps::ObservableSet>&, alps::Parameters const&) const {}

  template<typename WORKER>
  bool register_worker(std::string const& name) {
    bool isnew = (worker_creators_.find(name) == worker_creators_.end());
    worker_creators_[name] = worker_pointer_type(new parallel_worker_creator<WORKER>());
    return isnew;
  }

  bool unregister_worker(std::string const& name) {
    worker_map_type::iterator itr = worker_creators_.find(name);
    if (itr == worker_creators_.end()) return false;
    worker_creators_.erase(itr);
    return true;
  }

  static parallel_worker_factory* instance() {
    if (!instance_) instance_ = new parallel_worker_factory;
    return instance_;
  }

private:
  parallel_worker_factory() {}
  ~parallel_worker_factory() {}

  static parallel_worker_factory* instance_;

  typedef boost::shared_ptr<abstract_parallel_worker_creator> worker_pointer_type;
  typedef std::map<std::string, worker_pointer_type> worker_map_type;
  worker_map_type worker_creators_;
};
