/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004-2005 by Synge Todo <wistaria@comp-phys.org>
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
#include <alps/parameterlist.h>
#include <alps/lattice.h>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace looper;
#endif

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::ParameterList params(std::cin);

  for (int j = -1; j <= 1; ++j) {
    std::cout << "[Jxy = -1; Jz = " << j << "]\n";

    for (alps::ParameterList::const_iterator itr = params.begin();
         itr != params.end(); ++itr) {
      for (alps::Parameters::const_iterator p = itr->begin();
           p != itr->end(); ++p) {
        if (p->key() != "LATTICE_LIBRARY")
          std::cout << p->key() << " = " << p->value() << std::endl;
      }

      alps::graph_helper<> gh(*itr);
      looper::model_parameter model(gh.graph(), alps::half_integer<int>(0.5),
                                    -1.0, (double)j);

      std::cout << "model has "
                << (model.is_signed() ? "" : "no ")
                << "sign problem.\n";
      std::cout << "model is "
                << (model.is_frustrated() ? "" : "not ")
                << "classically frustrated.\n";
    }
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
