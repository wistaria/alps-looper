/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/temperature.h>
#include <alps/parameterlist.h>
#include <iostream>

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::ParameterList params(std::cin);

  for (alps::ParameterList::const_iterator p = params.begin();
       p != params.end(); ++p) {
    int steps = p->value_or_default("STEPS", 10);
    looper::temperature temp(*p);

    std::cout << "initial temperature = " << temp.initial() << std::endl
              << "final temperature   = " << temp.final() << std::endl
              << "annealing steps     = " << temp.annealing_steps()
              << std::endl;

    for (int i = 1; i <= steps; ++i)
      std::cout << i << ' ' << temp(i) << std::endl;
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
