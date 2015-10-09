/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/model_parameter.h>
#include <iostream>

template<typename G>
void output(const alps::graph_helper<G>& gh, const looper::model_parameter& m) {
  BOOST_FOREACH(typename alps::graph_helper<G>::site_descriptor s, gh.sites())
    std::cout << "site " << s << ": type = " << gh.site_type(s)
              << ": (S, C, Hx, Hz, D) = " << m.site(s, gh.graph()) << std::endl;
  BOOST_FOREACH(typename alps::graph_helper<G>::bond_descriptor b, gh.bonds())
    std::cout << "bond " << b << ": type = " << gh.bond_type(b)
              << ": (C, Jx, Jy, Jz) = " << m.bond(b, gh.graph()) << std::endl;
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
  looper::model_parameter mp(params, gh, mh);
  output(gh, mp);

  // construct for inhomogeneous model
  params["LOOPER_DEBUG[MODEL OUTPUT]"] = "cout";
  mp.set_parameters(params, gh.graph(), true, true, mh.model());

  // construct from parameters
  mp.set_parameters(gh.graph(), looper::site_parameter(1.5),
    looper::bond_parameter_xyz(0, -2, -2, -1));
  output(gh, mp);

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
