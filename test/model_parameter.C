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

#include <looper/model.h>
#include <boost/foreach.hpp>
#include <iostream>

std::ostream& operator<<(std::ostream& os, const looper::site_parameter& p) {
  os << "C = " << p.c << ", Hx = " << p.hx << ", Hz = " << p.hz;
  return os;
}

std::ostream& operator<<(std::ostream& os, const looper::bond_parameter& p) {
  os << "C = " << p.c << ", Jxy = " << p.jxy << ", Jz = " << p.jz;
  return os;
}

template<typename G>
void output(const alps::graph_helper<G>& gh, const looper::model_parameter& m) {
  BOOST_FOREACH(typename alps::graph_helper<G>::site_descriptor s, gh.sites())
    std::cout << "site " << s << ": type = " << gh.site_type(s)
              << ", S = " << m.site(s, gh.graph()).s << std::endl;
  BOOST_FOREACH(typename alps::graph_helper<G>::bond_descriptor b, gh.bonds())
    std::cout << "bond " << b << ": type = " << gh.bond_type(b)
              << ", " << m.bond(b, gh.graph()) << std::endl;
}

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::Parameters params;
  std::cin >> params;

  // get graph
  alps::graph_helper<> gh(params);

  // get model
  alps::model_helper<> mh(gh, params);

  // construct from model library
  looper::model_parameter m0(params, gh, mh);
  output(gh, m0);

  // construct for inhomogeneous model
  looper::model_parameter m1(params, gh.graph(), true, true, mh.model());
  output(gh, m1);

  // construct from parameters
  looper::model_parameter m2(gh.graph(), alps::half_integer<int>(1.5), -2, -1);
  output(gh, m2);

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
