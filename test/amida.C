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

/* $Id: amida.C 693 2004-03-16 15:48:04Z wistaria $ */

#include <looper/amida.h>
#include <looper/vector_helper.h>

#include <boost/random.hpp>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

const std::string dump_xdr = "amida.xdr";

const std::size_t n_series = 4;
const std::size_t n_nodes = 50;

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace looper;
#endif

template<class T>
void erase(looper::amida<T>& a, int n)
{
  looper::amida_node<T> * ptr =
    const_cast<looper::amida_node<T> *>(&(a.base()[n]));
  typename looper::amida<T>::iterator itr(ptr, ptr->series[0]);
  a.erase(itr);
}

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef boost::mt19937 rng_type;
  typedef boost::uniform_int<> uniform_int;
  rng_type rng;

  std::cout << "[[amida test]]\n";

  // amida array
  looper::amida<std::size_t> amida(n_series);

  // insert nodes randomly
  std::size_t n0, n1;
  for (std::size_t i = 0; i < n_nodes; i++) {
    n0 = uniform_int(0, n_series - 1)(rng);
    n1 = n0;
    while (n0 == n1) n1 = uniform_int(0, n_series - 1)(rng);
    looper::amida<std::size_t>::iterator r0 = amida.series(n0).first;
    looper::amida<std::size_t>::iterator r1 = amida.series(n1).first;
    looper::amida<std::size_t>::iterator r = 
      amida.insert_link_next(0, r0, r1).first;
    std::cout << "node connecting " << n0 << " and " << n1 
              << " inserted at " << index(amida, r) << std::endl;
  }

  erase(amida, 10);
  erase(amida, 16);
  erase(amida, 17);
  erase(amida, 12);

  std::cout << std::endl;
  std::cout << amida;

  std::cout << "\n[[ trace along series ]]\n";

  for (std::size_t i = 0; i < n_series; i++) {
    std::cout << i << '\t';
    looper::amida<std::size_t>::iterator a = amida.series(i).first;
    while (! a.at_top()) {
      std::cout << "(" << index(amida, a) << ',' << a.leg() << ")->";
      ++a;
    }
    std::cout << "(" << index(amida, a) << ',' << a.leg() << ")\n";
  }

  std::cout << "\n[[ trace along world line ]]\n";
  for (std::size_t i = 0; i < n_series; i++) {
    looper::amida<std::size_t>::iterator a = amida.series(i).first;
    std::cout << "[" << index(amida, a) << "," << a.series() << "]";
    ++a;
    while (! a.at_boundary()) {
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      a.jump();
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      ++a;
    }
    std::cout << "->[" << index(amida, a) << "," << a.series() << "]\n";
  }

  std::cout << "\n[[ trace along loop ]]\n";
  std::size_t d;

  for (std::size_t i = 0; i < n_series; i++) {
    looper::amida<std::size_t>::iterator a = amida.series(i).first;
    d = 0; // upward
    std::cout << "[" << index(amida, a) << "," << a.series() << "]";
    if (d == 0) ++a; else --a; // proceed(a, d);
    while (! a.at_boundary()) {
      d ^= 1;
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      a.jump();
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      if (d == 0) ++a; else --a; // proceed(a, d);
    }
    std::cout << "->[" << index(amida, a) << "," << a.series() << "]\n";
  }

  for (std::size_t i = 0; i < n_series; i++) {
    looper::amida<std::size_t>::iterator a = amida.series(i).second;
    d = 1; // downward
    std::cout << "[" << index(amida, a) << "," << a.series() << "]";
    if (d == 0) ++a; else --a; // proceed(a, d);
    while (! a.at_boundary()) {
      d ^= 1;
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      a.jump();
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      if (d == 0) ++a; else --a; // proceed(a, d);
    }
    std::cout << "->[" << index(amida, a) << "," << a.series() << "]\n";
  }

  std::cout << "[Save & load test]\n";

  {
    alps::OXDRFileDump dump(dump_xdr);
    std::cout << "Saving to " << dump_xdr << "...\n";
    dump << amida;
  }

  amida.clear();

  {
    alps::IXDRFileDump dump(dump_xdr);
    std::cout << "Loading from " << dump_xdr << "...\n";
    dump >> amida;
  }

  for (std::size_t i = 0; i < n_series; i++) {
    looper::amida<std::size_t>::iterator a = amida.series(i).second;
    d = 1; // downward
    std::cout << "[" << index(amida, a) << "," << a.series() << "]";
    if (d == 0) ++a; else --a; // proceed(a, d);
    while (! a.at_boundary()) {
      d ^= 1;
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      a.jump();
      std::cout << "->(" << index(amida, a) << "," << a.series() << ")";
      if (d == 0) ++a; else --a; // proceed(a, d);
    }
    std::cout << "->[" << index(amida, a) << "," << a.series() << "]\n";
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
