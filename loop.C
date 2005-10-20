/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef WITHOUT_SCHEDULER

#include <alps/osiris.h>
#include <alps/scheduler.h>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    return alps::scheduler::start(argc, argv, factory());
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

#else

#include <sstream>
#include <boost/timer.hpp>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  alps::Parameters params;
  switch (argc) {
  case 1:
    std::cin >> params;
    break;
  case 2:
    {
      boost::filesystem::ifstream is(argv[1]);
      is >> params;
      break;
    }
  default:
    boost::throw_exception(std::invalid_argument(
      "Usage: " + std::string(argv[0]) + " [paramfile]"));
  }

  std::cout << "[input parameters]\n" << params;

  qmc_worker_base<>* sim =
    factory().make_worker(alps::ProcessList(1), params, 0);
  boost::timer tm;

  bool thermalized = false;
  while (sim->work_done() < 1.0) {
    sim->dostep();
    if (!thermalized && sim->is_thermalized()) {
      const_cast<alps::ObservableSet&>(sim->get_measurements()).reset(true);
      thermalized = true;
    }
  }

  double t = tm.elapsed();
  std::cerr << "[speed]\nelapsed time = " << t << " sec ("
            << (sim.mcs()) / t << " MCS/sec)\n";

  std::cout << "[results]\n" << sim.get_measurements();

#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
    std::cerr << exc.what() << "\n";
    return -1;
  }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif
}

#endif
