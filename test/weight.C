/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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
#include <looper/weight.h>
#include <alps/parameterlist.h>
#include <iostream>

void output(const looper::site_parameter& p, const looper::site_weight& w)
{
  std::cout << "C = " << p.c()
	    << "Hx = " << p.hx()
	    << "Hz = " << p.hz()
	    << ": v[1] = " << w.v(1)
	    << ", v[2] = " << w.v(2)
	    << ", v[3] = " << w.v(3)
	    << ", offset = " << w.offset()
            << ", sign = " << weight.sign() << std::endl;
  looper::site_weight::check(p, w);
}

void output(const looper::bond_parameter& p, const looper::bond_weight& w)
{
  std::cout << "C = " << p.c()
	    << ", Jxy = " << p.jxy()
            << ", Jz = " << p.jz()
	    << ": v[1] = " << w.v(1)
	    << ", v[2] = " << w.v(2)
	    << ", v[3] = " << w.v(3)
	    << ", v[4] = " << w.v(4)
	    << ", offset = " << w.offset()
	    << ", sign = " << w.sign() << std::endl;
  looper::bond_weight::check(p, w);
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::ParameterList params;
  std::cin >> params;

  for (alps::ParameterList::iterator p = params.begin();
       p != params.end(); ++p) {
    looper::site_parameter site(0.5, (*p)["Hz"], (*p)["Hx"]);
    looper::bond_parameter bond((*p)["C"], (*p)["Jxy"], (*p)["Jz"]);

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
