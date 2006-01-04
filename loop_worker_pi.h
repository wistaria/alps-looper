/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOP_WORKER_PI_H
#define LOOP_WORKER_PI_H

#include "loop_worker.h"
#include <looper/cluster.h>
#include <looper/operator.h>
#include <looper/type.h>

class qmc_worker_pi : public qmc_worker_base
{
public:
  typedef looper::path_integral                          qmc_type;
  typedef qmc_worker_base                                super_type;

  typedef looper::local_operator<qmc_type, loop_graph_t> local_operator_t;
  typedef std::vector<local_operator_t>                  operator_string_t;
  typedef operator_string_t::iterator                    operator_iterator;

  typedef looper::union_find::node                       cluster_fragment_t;
  typedef looper::cluster_info                           cluster_info_t;

  qmc_worker_pi(alps::ProcessList const& w, alps::Parameters const& p, int n);
  virtual void dostep();
  void save(alps::ODump& dp) const;
  void load(alps::IDump& dp);

protected:
  template<typename BIPARTITE, typename FIELD, typename IMPROVE>
  void dostep_impl();

private:
  std::vector<int> spins;
  std::vector<local_operator_t> operators;

  // working area
  std::vector<int> spins_c;
  std::vector<local_operator_t> operators_p;
  std::vector<cluster_fragment_t> fragments;
  std::vector<int> current;
  std::vector<cluster_info_t> clusters;
  std::vector<improved_estimator_t::estimate::type> estimates;
};

#endif // LOOP_WORKER_PI_H
