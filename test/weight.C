/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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
#include <looper/weight_impl.h>
#ifdef HAVE_PARAPACK_13
# include <alps/math.hpp>
#else
# include <alps/numeric/round.hpp>
#endif
#include <boost/random.hpp>
#include <iostream>

#ifdef HAVE_PARAPACK_13
using alps::round;
#else
using alps::numeric::round;
#endif

void output(const looper::site_parameter& p, const looper::site_weight_helper& w)
{
  std::cout << "C = " << p.c
            << ", Hx = " << p.hx
            << " : v[0] = " << round<1>(w.v[0])
            << ", offset = " << w.offset
            << ", sign = " << w.sign << std::endl;
  w.check(p);
}

void output(const looper::bond_parameter_xxz& p, const looper::xxz_bond_weight_helper& w)
{
  std::cout << "C = " << p.c
            << ", Jxy = " << p.jxy
            << ", Jz = " << p.jz
            << " : v[0] = " << round<1>(w.v[0])
            << ", v[1] = " << round<1>(w.v[1])
            << ", v[2] = " << round<1>(w.v[2])
            << ", v[3] = " << round<1>(w.v[3])
            << ", offset = " << w.offset
            << ", sign = " << w.sign << std::endl;
  w.check(p);
}

void output(const looper::bond_parameter_xyz& p, const looper::xyz_bond_weight_helper& w)
{
  std::cout << "C = " << p.c
            << ", Jx = " << p.jx
            << ", Jy = " << p.jy
            << ", Jz = " << p.jz
            << " : v[0] = " << round<1>(w.v[0])
            << ", v[1] = " << round<1>(w.v[1])
            << ", v[2] = " << round<1>(w.v[2])
            << ", v[3] = " << round<1>(w.v[3])
            << ", v[4] = " << round<1>(w.v[4])
            << ", v[5] = " << round<1>(w.v[5])
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

  alps::Parameters p0;
  alps::Parameters p1;
  alps::Parameters p2;
  alps::Parameters p3;
  p1["FORCE_SCATTER"] = 0.1;
  p3["FORCE_SCATTER"] = 0.1;
  p2["SCATTERING_RATIO"] = 1.0;
  p3["SCATTERING_RATIO"] = 1.0;

  for (alps::ParameterList::const_iterator p = params.begin();
       p != params.end(); ++p) {
    looper::site_parameter site(0.5, 0, p->value_or_default("Hx",0), 0, 0);
    looper::bond_parameter_xxz bond(p->value_or_default("C",0),
      p->value_or_default("Jxy",0), p->value_or_default("Jz",0));

    std::cout << "site weight (standard): ";
    output(site, looper::site_weight_helper(site, p0));

    std::cout << "bond weight (standard): ";
    output(bond, looper::xxz_bond_weight_helper(bond, p0));

    std::cout << "bond weight (ergodic): ";
    output(bond, looper::xxz_bond_weight_helper(bond, p1));
  }

  std::cout << "[random check for xxz_bond_weight]\n";

  // random number generator
  boost::mt19937 eng(29833u);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    rng(eng, boost::uniform_real<>(-1,1));

  for (int i = 0; i < 10; ++i) {
    double sc = rng();
    double hx = rng();
    looper::site_parameter site(0.5, sc, hx, 0, 0);
    double bc = rng();
    double jxy = rng();
    double jz = rng();
    looper::bond_parameter_xxz bond(bc, jxy, jz);

    std::cout << "site weight (standard): ";
    output(site, looper::site_weight_helper(site, p0));

    std::cout << "bond weight (standard): ";
    output(bond, looper::xxz_bond_weight_helper(bond, p0));

    std::cout << "bond weight (ergodic): ";
    output(bond, looper::xxz_bond_weight_helper(bond, p1));
  }

  std::cout << "[random check for xyz_bond_weight]\n";

  for (int i = 0; i < 10; ++i) {
    double sc = rng();
    double hx = rng();
    looper::site_parameter site(0.5, sc, hx, 0, 0);
    double bc = rng();
    double jx = rng();
    double jy = std::abs(jx) * rng(); // |jy| should be smaller than |jx|
    double jz = rng();
    looper::bond_parameter_xyz bond(bc, jx, jy, jz);

    std::cout << "site weight (standard): ";
    output(site, looper::site_weight_helper(site, p0));

    std::cout << "bond weight (standard): ";
    output(bond, looper::xyz_bond_weight_helper(bond, p0));

    std::cout << "bond weight (ergodic): ";
    output(bond, looper::xyz_bond_weight_helper(bond, p1));

    std::cout << "bond weight (scattering): ";
    output(bond, looper::xyz_bond_weight_helper(bond, p2));

    std::cout << "bond weight (ergodic + scattering): ";
    output(bond, looper::xyz_bond_weight_helper(bond, p3));
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
