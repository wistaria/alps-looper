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

// $Id: fill_duration.C 554 2003-11-12 02:36:24Z wistaria $

#include <looper/fill_duration.h>
#include <alps/alea.h>
#include <boost/random.hpp>
#include <iostream>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  double r, tmax;
  int trial;
  std::cin >> r >> tmax >> trial;

  boost::mt19937 base_rng;
  boost::uniform_01<boost::mt19937> uniform_01(base_rng);
  for (int i = 0; i < trial; ++i) uniform_01();

  std::cout << "filling duration [0," << tmax << "]\n";
  std::vector<double> a;
  looper::fill_duration(uniform_01, a, r, tmax);
  for (std::vector<double>::iterator itr = a.begin(); itr != a.end(); ++itr) {
    std::cout << *itr << ' ';
  }
  std::cout << std::endl;
  
  std::cout << "filling duration " << trial << " times\n";
  alps::BasicSimpleObservable<double, alps::NoBinning<double> >
    n("average density");
  n.reset(true);
  for (int i = 0; i < trial; ++i) {
    looper::fill_duration(uniform_01, a, r, tmax);
    n << a.size() / tmax;
  }
  std::cout << n;
  std::cout << "exact density: " << r << std::endl;

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
