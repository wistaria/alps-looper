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

#include "loop_worker.h"

#include <alps/alea.h>
#include <alps/model.h>
#include <alps/scheduler.h>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <string>

namespace po = boost::program_options;

struct options
{
  boost::uint32_t seed;
  unsigned int dim;
  unsigned int lsize;
  double spin;
  double jxy;
  double jz;
  double temp;
  unsigned int step_t;
  unsigned int step_m;
  double fs;
};


int main(int argc, char *argv[])
{

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  // options
  options opts;
  po::options_description desc("Command line options");
  desc.add_options()
    ("seed,r", po::value<boost::uint32_t>(&opts.seed)->default_value(2837),
     "seed for random number generator")
    ("dim,d", po::value<unsigned int>(&opts.dim)->default_value(1),
     "spatial dimensions")
    ("size,l", po::value<unsigned int>(&opts.lsize)->default_value(16),
     "linear size of the system")
    ("spin,s", po::value<double>(&opts.spin)->default_value(0.5), "spin size S")
    ("jxy,x", po::value<double>(&opts.jxy)->default_value(1.0),
     "coupling Jxy (positive for antiferromagnetc)")
    ("jz,z", po::value<double>(&opts.jz)->default_value(1.0),
     "coupling Jz (positive for antiferromagnetic)")
    ("temp,t", po::value<double>(&opts.temp)->default_value(1.0), "temperature")
    ("therm,m", po::value<unsigned int>(&opts.step_t)->default_value(1024),
     "MCS for thermalization")
    ("measure,n", po::value<unsigned int>(&opts.step_m)->default_value(8192),
     "MCS for measurement")
    ("sse,e", "use SSE representation instead of path-integral")
    ("scatter,c", po::value<double>(&opts.fs)->default_value(0.0),
     "ratio of scattering graph (only for experts)")
    ("version,v", "print version and exit")
    ("help,h", "this help")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }  
  if (vm.count("version")) {
    looper::print_copyright(std::cout);
    alps::print_copyright(std::cout);
    return 0;
  }
  looper::print_copyright(std::cerr);
  alps::print_copyright(std::cerr);

  // simulation parameters
  alps::half_integer<int> spin(opts.spin);
  std::cout << "r:   seed for RNG           : " << opts.seed << std::endl
	    << "d:   spatial dimension      : " << opts.dim << std::endl
            << "l:   linear size            : " << opts.lsize << std::endl
            << "s:   spin size S            : " << spin << std::endl
            << "x:   coupling Jxy           : " << opts.jxy << std::endl
            << "z:   coupling Jz            : " << opts.jz << std::endl
            << "t:   temperature            : " << opts.temp << std::endl
            << "m:   MCS for thermalization : " << opts.step_t << std::endl
            << "n:   MCS for measurement    : " << opts.step_m << std::endl
            << "e:   representation         : "
	    << (vm.count("sse") == 0 ? "path integral" : "SSE")
            << std::endl << std::endl;

  // random number generator
  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    rng(boost::mt19937(opts.seed), boost::uniform_real<>());
  for (int i = 0; i < 19844; ++i) rng();

  // hypercubic lattice (real lattice)
  typedef looper::parity_graph_type graph_type;
  graph_type g;
  looper::generate_graph(g,
    looper::hypercubic_graph_generator<>(opts.dim, opts.lsize));

  // model & inverse temperature
  typedef looper::model_parameter<> model_type;
  model_type model(-opts.jxy, -opts.jz, spin, g, false);
  double beta = 1./opts.temp;

  // measurements
  alps::ObservableSet measurements;

  if (vm.count("sse") == 0) {

    // path-integral representation
    qmc_worker<looper::path_integral<graph_type,
      model_type> > worker(g, model, beta, opts.fs, measurements);

    for (int mcs = 0; mcs < opts.step_t + opts.step_m; ++mcs) {
      if (mcs == opts.step_t) measurements.reset(true);
      worker.step(rng, measurements);
    }
    accumulate(measurements, measurements);
    looper::print_all(std::cout, measurements);

  } else {

    // SSE representation
    qmc_worker<looper::sse<graph_type,
      model_type> > worker(g, model, beta, opts.fs, measurements);

    for (int mcs = 0; mcs < opts.step_t + opts.step_m; ++mcs) {
      if (mcs == opts.step_t) measurements.reset(true);
      worker.step(rng, measurements);
    }
    accumulate(measurements, measurements);
    looper::print_all(std::cout, measurements);

  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (const std::exception& excp) {
  std::cerr << excp.what() << std::endl;
  std::exit(-1); }
catch (...) {
  std::cerr << "Unknown exception occurred!" << std::endl;
  std::exit(-1); }
#endif
}
