/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2003-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include "loop_config.h"
#include <looper/evaluate.h>
#include <alps/alea.h>
#include <alps/scheduler.h>
#include <boost/filesystem/operations.hpp>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef looper::evaluator<loop_config::estimator_t> evaluator_t;

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " file1 [file2 [...]]\n";
    std::exit(-1);
  }

  alps::scheduler::SimpleMCFactory<evaluator_t> evaluator_factory;
  alps::scheduler::init(evaluator_factory);
  for (int i = 1; i < argc; ++i) {
    boost::filesystem::path p = complete(boost::filesystem::path(argv[i]));
    alps::ProcessList nowhere;
    alps::scheduler::MCSimulation sim(nowhere, p);
    if (sim.runs.size()) {
      dynamic_cast<evaluator_t*>(sim.runs[0])->evaluate(sim);
      sim.checkpoint(p);
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif
}
