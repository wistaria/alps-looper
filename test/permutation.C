/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/permutation.h>
#include <boost/random.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

const unsigned int n = 10;
const unsigned int trial1 = 20;
const unsigned int trial2 = 5;
const boost::uint32_t seed = 23094;

int main() {
  typedef boost::mt19937 rng_type;
  typedef std::vector<unsigned int> vector_type;

  rng_type rng(seed);
  boost::variate_generator<rng_type&, boost::uniform_real<> >
    random(rng, boost::uniform_real<>(0, 1));

  std::cout << "[[random permutaion test 1]]\n";

  std::cout << "generating " << trial1 << " random permutations of "
       << n << " integers [0..." << n - 1 << "]\n";

  for (unsigned int i = 0; i < trial1; ++i) {
    vector_type result(n);
    for (unsigned int j = 0; j < n; ++j) result[j] = j;
    looper::random_shuffle(result.begin(), result.end(), random);
    for (unsigned int j = 0; j < n; ++j) {
      std::cout << result[j] << '\t';
    }
    std::cout << std::endl;
  }

  std::cout << "\n[[random permutaion test 2 (with restrictions)]]\n";

  std::cout << "generating " << trial2 << " partitioned random permutations\n";

  for (unsigned int t = 0; t < trial2; t++) {
    vector_type config0(n);
    for (unsigned int i = 0; i < n; i++) {
      config0[i] = boost::uniform_smallint<>(0, 1)(rng);
    }
    vector_type config1(config0);
    looper::random_shuffle(config1.begin(), config1.end(), random);

    std::cout << "conf0";
    for (unsigned int i = 0; i < n; i++) std::cout << '\t' << config0[i];
    std::cout << std::endl;

    std::cout << "conf1";
    for (unsigned int i = 0; i < n; i++) std::cout << '\t' << config1[i];
    std::cout << std::endl;

    vector_type result(n);
    for (unsigned int j = 0; j < n; ++j) result[j] = j;

    // generate partitioned permutation
    looper::partitioned_random_shuffle(result.begin(), result.end(),
                                       config0.begin(), config1.begin(),
                                       random);

    std::cout << "perm";
    for (unsigned int i = 0; i < n; i++) std::cout << '\t' << result[i];
    std::cout << std::endl;

    std::cout << "check";
    int chk = 0;
    for (unsigned int i = 0; i < n; i++) {
      std::cout << '\t' << (config0[i] ^ config1[result[i]]);
      chk += (config0[i] ^ config1[result[i]]);
    }
    std::cout << std::endl;
    if (chk > 0) { std::cout << "Error occurs!  Stop.\n"; std::exit(-1); }

    std::cout << std::endl;
  }
}
