/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004 by Synge Todo <wistaria@comp-phys.org>
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

#include <cmath>
#include <boost/random.hpp>
#include <algorithm>
#include <iostream>
#include <vector>

// Intel C++ (icc) 8.0 fails in this test when -xW vectorized option
// is specified.  -- ST 2004.04.01

template<class RNG>
void generate(RNG& uniform_01)
{
  typedef RNG rng_type;
  std::vector<int> array(20);
  std::generate(array.begin(), array.end(),
                boost::variate_generator<rng_type&, boost::uniform_smallint<> >
                (uniform_01, boost::uniform_smallint<>(0, 1)));

  std::vector<int>::const_iterator itr_end = array.end();
  for (std::vector<int>::const_iterator itr = array.begin(); itr != itr_end;
       ++itr) std::cout << *itr << ' ';
  std::cout << std::endl;
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    uniform_01(boost::mt19937(4357), boost::uniform_real<>());
  generate(uniform_01);

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
