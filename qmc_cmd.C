/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

#include "qmc_impl.h"
#include <looper/exact_diag.h>

#include <alps/alea.h>
#include <boost/random.hpp>
#include <iostream>
#include <string>

struct options {
  uint32_t seed;                // -r seed
  uint32_t dim;                 // -d spatial dimension
  uint32_t lsize;               // -l linear size of system
  alps::half_integer<int> spin; // -s spin size S
  double Jxy;                   // -x coupling Jxy (positive for ferromagnetic)
  double Jz;                    // -z coupling Jz (positive for ferromagnetic)
  double temp;                  // -t temperature
  uint32_t step_t;              // -m MCS for thermalization
  uint32_t step_m;              // -n MCS for measurement
  std::string representation;   // -e use SSE instead of path integral
                                // -g use exact diagonalization

  options(int argc, char *argv[]) : 
    // default options
    seed(2837), dim(1), lsize(16), spin(0.5), Jxy(-1.), Jz(-1.), temp(1.),
    step_t(1024), step_m(8192), representation("path integral")
  {
    parse(argc, argv);
  }
  
  void usage(int status, std::ostream& os = std::cerr) const {
    os << "[command line options]\n\n"
       << "  -r int     seed of random number generator\n"
       << "  -d int     spatial dimension\n"
       << "  -l int     linear size of system\n"
       << "  -s double  spin size S\n"
       << "  -x double  coupling Jxy\n"
       << "  -z double  coupling Jz\n"
       << "  -t double  temperature\n"
       << "  -m int     MCS for thermalization\n"
       << "  -n int     MCS for measurement\n"
       << "  -e         use SSE representation instead of path-integral one\n"
#ifdef HAVE_LAPACK
       << "  -g         use exact diagonalization instead of QMC\n"
#endif // HAVE_LAPACK
       << "  -h         this help\n\n";
    if (status) {
      boost::throw_exception(std::invalid_argument("Invalid command line option(s)"));
    } else {
      std::exit(0);
    }
  }

  void parse(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'r' :
	  if (i + 1 == argc) usage(1);
          seed = std::atoi(argv[++i]);
          break;
        case 'd' :
	  if (i + 1 == argc) usage(1);
          dim = std::atoi(argv[++i]);
          break;
        case 'l' :
	  if (i + 1 == argc) usage(1);
          lsize = std::atoi(argv[++i]);
          break;
        case 's' :
	  if (i + 1 == argc) usage(1);
          spin = alps::half_integer<int>(std::atof(argv[++i]));
          break;
        case 'x' :
	  if (i + 1 == argc) usage(1);
          Jxy = std::atof(argv[++i]);
          break;
        case 'z' :
	  if (i + 1 == argc) usage(1);
          Jz = std::atof(argv[++i]);
          break;
        case 't' :
	  if (i + 1 == argc) usage(1);
          temp = std::atof(argv[++i]);
          break;
        case 'm' :
	  if (i + 1 == argc) usage(1);
          step_t = std::atoi(argv[++i]);
          break;
        case 'n' :
	  if (i + 1 == argc) usage(1);
          step_m = std::atoi(argv[++i]);
          break;
        case 'e' :
	  representation = "SSE";
          break;
#ifdef HAVE_LAPACK
        case 'g' :
	  representation = "exact diagonalization";
          break;
#endif // HAVE_LAPACK
	case 'h' :
	  usage(0);
	  break;
	default :
	  usage(1);
	  break;
	}
	break;
	
      default :
	usage(1);
	break;
      }
    }
  }

  void output(std::ostream& os = std::cout) const
  {
    os << seed << ' '
       << dim << ' '
       << lsize << ' '
       << spin << ' '
       << Jxy << ' '
       << Jz << ' '
       << temp << ' '
       << step_t << ' '
       << step_m << ' ';
    if (representation == "path integral") {
      os << "PI";
    } else if (representation == "SSE") {
      os << "SSE";
#ifdef HAVE_LAPACK
    } else {
      os << "ED";
#endif // HAVE_LAPACK
    }
  }
};


int main(int argc, char *argv[])
{

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  options opts(argc, argv);
  looper::print_copyright(std::cout);
  alps::print_copyright(std::cout);
  std::cout <<
    "qmc_cmd: a command-line QMC program for simple hypercubic lattices\n\n";

  // simulation parameters
  std::cout << "r:   seed for RNG           : " << opts.seed << std::endl
	    << "d:   spatial dimension      : " << opts.dim << std::endl
	    << "l:   linear size            : " << opts.lsize << std::endl
	    << "s:   spin size S            : " << opts.spin << std::endl
	    << "x:   coupling Jxy           : " << opts.Jxy << std::endl
	    << "z:   doupling Jz            : " << opts.Jz << std::endl
	    << "t:   temperature            : " << opts.temp << std::endl
	    << "m:   MCS for thermalization : " << opts.step_t << std::endl
	    << "n:   MCS for measurement    : " << opts.step_m << std::endl
#ifdef HAVE_LAPACK
	    << "e,g: representation         : " << opts.representation
#else
	    << "e:   representation         : " << opts.representation
#endif // HAVE_LAPACK
	    << std::endl << std::endl;

  // random number generator
  boost::mt19937 base_rng;
  boost::uniform_01<boost::mt19937> rng(base_rng);
  rng.base().seed(boost::mt19937::result_type(opts.seed));
  for (int i = 0; i < 19844; ++i) rng();

  // hypercubic lattice (real lattice)
  typedef looper::parity_graph_type graph_type;
  graph_type g;
  looper::generate_graph(
    looper::hypercubic_graph_generator<>(opts.dim, opts.lsize), g);

  // model & inverse temperature
  typedef looper::xxz_model model_type;
  model_type model(opts.Jxy, opts.Jz, opts.spin, g);
  double beta = 1./opts.temp;

  // measurements
  alps::ObservableSet measurements;

  if (opts.representation == "path integral" ||
      opts.representation == "SSE") {
    
    // random number generator
    boost::mt19937 base_rng;
    boost::uniform_01<boost::mt19937> rng(base_rng);
    rng.base().seed(boost::mt19937::result_type(opts.seed));
    for (int i = 0; i < 19844; ++i) rng();
    
    if (opts.representation == "path integral") {
      // path-integral representation
      qmc_worker<looper::path_integral<looper::virtual_graph<graph_type>,
	model_type> > worker(g, model, beta, measurements);
      
      for (int mcs = 0; mcs < opts.step_t + opts.step_m; ++mcs) {
	if (mcs == opts.step_t) measurements.reset(true);
	worker.step(rng, measurements);
      }
      worker.accumulate(measurements);

      std::cout << measurements << std::endl; 
      opts.output(); std::cout << ' ';
      worker.output_results(std::cout, measurements);
      std::cout << std::endl;
    } else {
      // SSE representation
      qmc_worker<looper::sse<looper::virtual_graph<graph_type>,
	model_type> > worker(g, model, beta, measurements);
      
      for (int mcs = 0; mcs < opts.step_t + opts.step_m; ++mcs) {
	if (mcs == opts.step_t) measurements.reset(true);
	worker.step(rng, measurements);
      }
      worker.accumulate(measurements);
      
      std::cout << measurements << std::endl; 
      opts.output(); std::cout << ' ';
      worker.output_results(std::cout, measurements);
      std::cout << std::endl;
    }
#ifdef HAVE_LAPACK
  } else {
    // exact diagonalization
    measurements.reset(true);
    ed_worker<looper::exact_diagonalization<graph_type, model_type> >
      worker(g, model, beta, measurements);
    worker.step(rng, measurements);
    worker.accumulate(measurements);

    std::cout << measurements << std::endl; 
    opts.output(); std::cout << ' ';
    worker.output_results(std::cout, measurements);
    std::cout << std::endl;
#endif // HAVE_LAPACK
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
