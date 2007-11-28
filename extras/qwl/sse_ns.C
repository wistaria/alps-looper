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

#include "sse.h"
#include <parapack/scheduler.h>
#include <boost/timer.hpp>
#include <time.h>

typedef loop_worker worker_type;
typedef alps::parapack::default_evaluator evaluator_type;

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  alps::ParameterList parameterlist;
  switch (argc) {
  case 1:
    std::cin >> parameterlist;
    break;
  case 2:
    {
      boost::filesystem::ifstream is(argv[1]);
      is >> parameterlist;
      break;
    }
  default:
    boost::throw_exception(std::invalid_argument(
      "Usage: " + std::string(argv[0]) + " [paramfile]"));
  }

  BOOST_FOREACH(alps::Parameters& p, parameterlist) {
    boost::timer tm;
    if (!p.defined("SEED")) p["SEED"] = static_cast<unsigned int>(time(0));
    if (!p.defined("WORKER_SEED")) p["WORKER_SEED"] = p["SEED"];
    if (!p.defined("DISORDER_SEED")) p["DISORDER_SEED"] = p["WORKER_SEED"];
    std::cout << "[input parameters]\n" << p;
    std::vector<alps::ObservableSet> obs;
    alps::parapack::abstract_worker* worker = new worker_type(p, obs);
    bool thermalized = worker->is_thermalized();
    while (worker->progress() < 1.0) {
      worker->run(obs);
      if (!thermalized && worker->is_thermalized()) {
        BOOST_FOREACH(alps::ObservableSet& s, obs) s.reset(true);
        thermalized = true;
      }
    }
    delete worker;

    BOOST_FOREACH(alps::ObservableSet& s, obs) evaluator_type::evaluate_observable(s, p);
    std::cerr << "[speed]\nelapsed time = " << tm.elapsed() << " sec\n";
    std::cout << "[results]\n";
    for (int i = 0; i < obs.size(); ++i)
      std::cout << "[[parition " << i << "]]\n" << obs[i];
  }

#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
    alps::comm_exit(true);
    return -1;
  }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}
