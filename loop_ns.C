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

#include "loop_factory.h"
#include <alps/osiris/comm.h>
#include <alps/scheduler.h>
#include <boost/timer.hpp>
#include <time.h>

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

  for (alps::ParameterList::iterator p = parameterlist.begin();
       p != parameterlist.end(); ++p) {
    boost::timer tm;
    if (!p->defined("SEED")) (*p)["SEED"] = static_cast<unsigned int>(time(0));
    std::cout << "[input parameters]\n" << *p;
    alps::scheduler::MCRun* worker =
      loop_factory::instance()->make_worker(alps::ProcessList(1), *p, 0);
    bool thermalized = worker->is_thermalized();
    while (worker->work_done() < 1.0) {
      worker->dostep();
      if (!thermalized && worker->is_thermalized()) {
        const_cast<alps::ObservableSet&>(worker->get_measurements()).reset(1);
        thermalized = true;
      }
    }
    alps::ObservableSet m = worker->get_measurements();
    delete worker;

    looper::abstract_evaluator* evaluator = loop_factory::instance()->make_evaluator(*p);
    evaluator->pre_evaluate(m, *p, m);
    evaluator->evaluate(m, *p, m);
    delete evaluator;

    std::cerr << "[speed]\nelapsed time = " << tm.elapsed() << " sec\n";
    std::cout << "[results]\n" << m;
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
