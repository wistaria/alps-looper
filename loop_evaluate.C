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
#include <boost/program_options.hpp>

class options
{
public:
  typedef std::vector<std::string>::iterator file_iterator;
  typedef std::vector<std::string>::const_iterator const_file_iterator;

  options() {}
  options(int argc, char** argv) { parse(argc, argv); }

  bool parse(int argc, char** argv)
  {
    namespace po = boost::program_options;

    po::options_description desc("options");
    desc.add_options()
      ("help", "produce help message")
      ("t_min", po::value<double>(), "set minimum temperature (T_MIN)")
      ("t_max", po::value<double>(), "set maximum temperature (T_MAX)")
      ("t_delta", po::value<double>(),
       "set step width of temperature (T_DELTA)")
      ("input-file", po::value<std::vector<std::string> >(&files_),
       "input file");

    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(desc).positional(p).run(),
          vm);
    notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return false;
    }

    if (vm.count("t_min")) params_["T_MIN"] = vm["t_min"].as<double>();
    if (vm.count("t_max")) params_["T_MAX"] = vm["t_min"].as<double>();
    if (vm.count("t_delta")) params_["T_DELTA"] = vm["t_delta"].as<double>();

    return true;
  }

  alps::Parameters const& parameters() const { return params_; }
  unsigned int num_files() const { return files_.size(); }
  boost::filesystem::path file(unsigned int i) const
  { return complete(boost::filesystem::path(files_[i])); }

private:
  alps::Parameters params_;
  std::vector<std::string> files_;
};

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  //  typedef looper::evaluator<loop_config::estimator_t> evaluator_t;

  options opts;
  if (!opts.parse(argc, argv)) std::exit(-1);

//   alps::scheduler::SimpleMCFactory<evaluator_t> evaluator_factory;
//   alps::scheduler::init(evaluator_factory);
//   for (int i = 0; i < opts.num_files(); ++i) {
//     boost::filesystem::path f = opts.file(i);
//     alps::ProcessList nowhere;
//     alps::scheduler::MCSimulation sim(nowhere, f);
//     if (sim.runs.size()) {
//       dynamic_cast<evaluator_t*>(sim.runs[0])->
//         evaluate(sim, opts.parameters(), f);
//       sim.checkpoint(p);
//     }
//   }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif
}
