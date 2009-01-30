/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/scheduler.h>
#include <alps/parapack/parallel.h>
#include <boost/timer.hpp>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  alps::ParameterList parameterlist;
  if (world.rank() == 0) {
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
      // default:
    }
    for (int p = 1; p < world.size(); ++p) world.send(p, 0, parameterlist);
  } else {
    world.recv(0, 0, parameterlist);
  }

  for (alps::ParameterList::iterator p = parameterlist.begin();
       p != parameterlist.end(); ++p) {
    world.barrier();
    boost::timer tm;
    if (!p->defined("SEED")) (*p)["SEED"] = static_cast<unsigned int>(time(0));
    (*p)["WORKER_SEED"] = static_cast<unsigned int>((*p)["SEED"]) ^ world.rank();
    (*p)["DISORDER_SEED"] = (*p)["WORKER_SEED"];
    if (world.rank() == 0) std::cout << "[input parameters]\n" << *p;
    alps::ObservableSet m;
    boost::shared_ptr<alps::parapack::abstract_worker> worker
      = alps::parapack::parallel_worker_factory::make_worker(world, *p);
    worker->init_observables(*p, m);
    bool thermalized = worker->is_thermalized();
    while (worker->progress() < 1.0) {
      worker->run(m);
      if (!thermalized && worker->is_thermalized()) {
        m.reset(true);
        thermalized = true;
      }
    }

    if (world.rank() == 0) {
      boost::shared_ptr<alps::parapack::abstract_evaluator> evaluator
        = alps::parapack::evaluator_factory::make_evaluator(*p);
      evaluator->evaluate(m);
      std::cerr << "[speed]\nelapsed time = " << tm.elapsed() << " sec\n";
      std::cout << "[results]\n" << m;
    }
  }

  world.barrier();

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
