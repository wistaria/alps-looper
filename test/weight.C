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

#include <looper/util.h>
#include <looper/xxz.h>
#include <looper/weight.h>
#include <alps/parameterlist.h>
#include <iostream>

template<typename W>
void output(const alps::Parameters& param, const W& weight)
{
  using looper::range_01;

  assert(weight.p_freeze_para() == weight.p_freeze(0,0) &&
	 weight.p_freeze_para() == weight.p_freeze(1,1));
  assert(weight.p_freeze_anti() == weight.p_freeze(0,1) &&
	 weight.p_freeze_anti() == weight.p_freeze(1,0));
  assert(weight.p_accept_para() == weight.p_accept(0,0) &&
	 weight.p_accept_para() == weight.p_accept(1,1));
  assert(weight.p_accept_anti() == weight.p_accept(0,1) &&
	 weight.p_accept_anti() == weight.p_accept(1,0));
  assert(weight.sign() == 1 || weight.sign() == -1);

  std::cout << "Jxy = " << param["Jxy"]
	    << ", Jz = " << param["Jz"]
	    << " : r = " << weight.density();
  if (weight.density() > 0) {
    std::cout << ", parallel: accept = " << range_01(weight.p_accept_para());
    if (weight.p_accept_para() > 0)
      std::cout << ", freeze = " << range_01(weight.p_freeze_para());
    std::cout << ", antiparallel: accept = "
	      << range_01(weight.p_accept_anti());
    if (weight.p_accept_anti() > 0)
      std::cout << ", freeze = " << range_01(weight.p_freeze_anti());
    std::cout << ", P_r = " << range_01(weight.p_reflect())
	      << ", sign = " << weight.sign();
  }
  std::cout << std::endl;
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
    looper::xxz_parameter xxz(0, (*p)["Jxy"], (*p)["Jz"]);
    looper::default_weight w(xxz);
    output(*p, w);
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
