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

// $Id: random_choice.C 554 2003-11-12 02:36:24Z wistaria $

// Define the following macro if you want the original initialization
// routine of O(N^2).

/* #define USE_INIT_WALKER1977 */

#include <looper/random_choice.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <vector>

static const unsigned int n = 8;
static const unsigned int samples = 100000;

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  std::cout << "number of bins = " << n << std::endl;
  std::cout << "number of samples = " << samples << std::endl;

  typedef boost::mt19937 base_rng_type;
  typedef boost::uniform_01<base_rng_type> rng_type;

  std::vector<double> weights(n);

  // random number generator
  base_rng_type base_rng(29411);
  rng_type rng(base_rng);

  // generate weights
  std::generate(weights.begin(), weights.end(), rng);
  double tw = 0;
  for (unsigned int i = 0; i < n; ++i) tw += weights[i];

  // random_choice
  looper::random_choice<> rc(weights);

  // check
  if (rc.check(weights)) {
    std::cout << "check succeeded\n";
  } else {
    std::cout << "check failed\n";
    std::exit(-1);
  }    

  std::vector<double> accum(n, 0);
  for (unsigned int t = 0; t < samples; ++t) ++accum[rc(rng)];

  std::cout << "bin\tweight\t\tresult\t\tdiff\t\tsigma\t\tdiff/sigma\n";
  for (unsigned int i = 0; i < n; ++i) {
    double diff = std::abs((weights[i] / tw) - (accum[i] / samples));
    double sigma = std::sqrt(accum[i]) / samples;
    std::cout << i << "\t" << (weights[i] / tw) << "    \t" << (accum[i] / samples)
	      << "    \t" << diff << "    \t" << sigma << "    \t" << (diff / sigma)
	      << std::endl;
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
