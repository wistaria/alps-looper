/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/lapack.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/random.hpp>
#include <cstdlib>
#include <iostream>

struct Options {
  uint32_t n;        // -n dimension of matrix
  uint32_t seed;     // -r seed
  
  Options(int argc, char *argv[]) : 
    // default options
    n(128), seed(2837)
  { parse(argc, argv); }
  
  void usage(int status, std::ostream& os = std::cerr) const {
    os << "[command line options]\n\n"
       << "  -n int     dimension of matrix\n"
       << "  -r int     seed for random number generator\n"
       << "  -h         this help\n\n";
    if (status) {
      os << "Invalid command line option(s)\n";
      std::exit(status);
    } else {
      std::exit(0);
    }
  }

  void parse(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'n' :
          if (i + 1 == argc) usage(1);
          n = std::atoi(argv[++i]);
          break;
        case 'r' :
          if (i + 1 == argc) usage(1);
          seed = std::atoi(argv[++i]);
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
};


int main(int argc, char ** argv)
{
  std::cout << "Diagonalization of random real symmetric matrix\n";
  std::cerr << "[starting program: "
            << boost::posix_time::microsec_clock::local_time() << "]\n";

  Options opts(argc, argv);
  std::cout << "n: dimension of matrix : " << opts.n << std::endl;
  std::cout << "r: seed for RNG        : " << opts.seed << std::endl;

  std::cout << "initialization... ";

  boost::mt19937 rng(uint32_t(opts.seed));
  double norm = 1.0/rng.max();
  for (int i = 0; i < 1000; ++i) rng();

  const int n = opts.n;
  boost::numeric::ublas::vector<double> vec(n);
  boost::numeric::ublas::matrix<double> mat(n, n);
  for (int i = 0; i < n; ++i)
    for (int j = i; j < n; ++j)
      mat(i, j) = norm * double(rng()) - 0.5;
  std::cout << "done\n";

  boost::posix_time::ptime t0 =
    boost::posix_time::microsec_clock::local_time();
  std::cerr << "[starting diagonalization: " << t0 << "]\n";

  std::cout << "diagonalization... " << std::flush;
  looper::diagonalize(mat, vec);
  std::cout << "done\n";

  boost::posix_time::ptime t1 =
    boost::posix_time::microsec_clock::local_time();
  std::cerr << "[finished diagonalization: " << t1 << "]\n";
  std::cerr << "[elapsed time: " << t1 - t0 << "]\n";

  std::cout << "selected eigenvalues:\n";
  bool print_last = false;
  for (int i = 0; i < n; i += std::max(n/10, 1)) {
    std::cout << "    e(" << i << ")\t= " << vec(i) << std::endl;
    if (i == n - 1) print_last = true;
  }
  if (!print_last)
    std::cout << "    e(" << n-1 << ")\t= " << vec(n-1) << std::endl;
}
