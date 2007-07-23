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
#include <looper/weight.h>
#include <alps/math.hpp>
#include <boost/random.hpp>
#include <iostream>

void output(const looper::site_parameter& p, const looper::site_weight& w)
{
  std::cout << "C = " << p.c
            << ", Hx = " << p.hx
            << " : v[0] = " << alps::round<1>(w.v[0])
            << ", offset = " << w.offset
            << ", sign = " << w.sign << std::endl;
  w.check(p);
}

void output(const looper::bond_parameter_xxz& p, const looper::bond_weight& w)
{
  std::cout << "C = " << p.c
            << ", Jxy = " << p.jxy
            << ", Jz = " << p.jz
            << " : v[0] = " << alps::round<1>(w.v[0])
            << ", v[1] = " << alps::round<1>(w.v[1])
            << ", v[2] = " << alps::round<1>(w.v[2])
            << ", v[3] = " << alps::round<1>(w.v[3])
            << ", offset = " << w.offset
            << ", sign = " << w.sign << std::endl;
  w.check(p);
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[standard input]\n";

  alps::ParameterList params;
  std::cin >> params;

  for (alps::ParameterList::const_iterator p = params.begin();
       p != params.end(); ++p) {
    looper::site_parameter site(0.5, 0, p->value_or_default("Hx",0), 0, 0);
    looper::bond_parameter_xxz bond(p->value_or_default("C",0),
      p->value_or_default("Jxy",0), p->value_or_default("Jz",0));

    std::cout << "site weight (standard): ";
    output(site, looper::site_weight(site));

    std::cout << "bond weight (standard): ";
    output(bond, looper::bond_weight(bond));

    std::cout << "bond weight (ergodic): ";
    output(bond, looper::bond_weight(bond, 0.1));
  }

  std::cout << "[random check]\n";

  // random number generator
  boost::mt19937 eng(29833u);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    rng(eng, boost::uniform_real<>(-1,1));

  for (int i = 0; i < 10; ++i) {
    // NOTE: gnu4 does not compile correctly the following code
    // looper::site_parameter site(0.5, rng(), rng(), rng());
    double sc = rng();
    double hx = rng();
    looper::site_parameter site(0.5, sc, hx, 0, 0);
    double bc = rng();
    double jxy = rng();
    double jz = rng();
    looper::bond_parameter_xxz bond(bc, jxy, jz);

    std::cout << "site weight (standard): ";
    output(site, looper::site_weight(site));

    std::cout << "bond weight (standard): ";
    output(bond, looper::bond_weight(bond));

    std::cout << "bond weight (ergodic): ";
    output(bond, looper::bond_weight(bond, 0.1));
  }

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
