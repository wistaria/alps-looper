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

/* $Id: percolation_cmd.C 693 2004-03-16 15:48:04Z wistaria $ */

// percolation_cmd - a command-line percolation program for simple
// hypercubic lattices

#include "percolation_impl.h"
#include <looper/graph.h>

#include <alps/alea.h>
#include <boost/random.hpp>
#include <iostream>

struct Options {
  uint32_t seed;         // -r seed
  bool     bond;         // -b whether bond percolation or not
  uint32_t dim;          // -d spatial dimension
  uint32_t lsize;        // -l linear size of system
  uint32_t nm;           // -m number of samples
  double   p;            // -p occupation probability
  
  Options(int argc, char *argv[]) : 
    // default options
    seed(2837), bond(false), dim(2), lsize(32), nm(128), p(-1)
  {
    parse(argc, argv);
  }
  
  void usage(int status, std::ostream& os = std::cerr) const {
    os << "[command line options]\n\n"
       << "  -p double  concentration (REQUIED)\n"
       << "  -r int     seed of random number generator (2837)\n"
       << "  -b         bond percolation (false)\n"
       << "  -d int     spatial dimension (2)\n"
       << "  -l int     linear size of system (32)\n"
       << "  -m int     number of samples (128)\n"
       << "  -h         this help\n\n";
    if (status) {
      boost::throw_exception(std::invalid_argument("Invalid command line option(s)"));
    } else {
      std::exit(0);
    }
  }

  void parse(int argc, char *argv[]) {
    bool found_p = false;
    for (int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'r' :
          if (i + 1 == argc) usage(1);
          seed = std::atoi(argv[++i]);
          break;
        case 'b' :
          bond = true;
          break;
        case 'd' :
          if (i + 1 == argc) usage(1);
          dim = std::atoi(argv[++i]);
          break;
        case 'l' :
          if (i + 1 == argc) usage(1);
          lsize = std::atoi(argv[++i]);
          break;
        case 'm' :
          if (i + 1 == argc) usage(1);
          nm = std::atoi(argv[++i]);
          break;
        case 'p' :
          if (i + 1 == argc) usage(1);
          p = std::atof(argv[++i]);
          found_p = true;
          break;
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
    if (!found_p) usage(1);
  }
};


int main(int argc, char *argv[]) {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  looper::print_copyright();
  std::cout << "percolation_cmd - a command-line percolation program for simple hypercubic lattices\n\n";

  // simulation parameters
  Options opts(argc, argv);
  std::cout << "r: seed for RNG           : " << opts.seed << std::endl;
  std::cout << "b: type                   : " << (opts.bond ? "bond" : "site") << " percolation\n";
  std::cout << "d: spatial dimension      : " << opts.dim << std::endl;
  std::cout << "l: linear size            : " << opts.lsize << std::endl;
  std::cout << "m: number of samples      : " << opts.nm << std::endl;
  std::cout << "p: concentration          : " << opts.p << std::endl
            << std::endl;

  // random number generator
  boost::mt19937 base_rng;
  boost::uniform_01<boost::mt19937> rng(base_rng);
  rng.base().seed(boost::mt19937::result_type(opts.seed));
  for (int i = 0; i < 19844; ++i) rng();

  // hypercubic lattice
  looper::graph_type g;
  looper::generate_graph(
    looper::hypercubic_graph_generator<>(opts.dim, opts.lsize), g);

  // measurements
  alps::ObservableSet measurements;

  percolation::worker_base wb(g, opts.p, opts.bond, measurements);
  for (int s = 0; s < opts.nm; s++) wb.step(g, rng, measurements);
  
  // output results
  std::cout << measurements;
  std::cout << std::endl
            << opts.seed << ' '
            << (opts.bond ? "bond" : "site") << ' '
            << opts.dim << ' '
            << opts.lsize << ' '
            << opts.nm << ' '
            << opts.p << ' ';
  wb.output_results(std::cout, measurements);
  std::cout << std::endl;

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
