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
#include <iostream>

template<typename G, typename SITE_P, typename BOND_P>
void output(const alps::graph_helper<G>& gh,
            const looper::model_parameter<SITE_P, BOND_P>& m)
{
  typedef G graph_type;
  typedef typename alps::graph_helper<G>::vertex_iterator vertex_iterator;
  typedef typename alps::graph_helper<G>::edge_iterator edge_iterator;

  // site parameters
  vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = gh.vertices(); vi != vi_end; ++vi) {
    std::cout << "site " << *vi << ": type = " << gh.vertex_type(*vi)
              << ", S = " << m.site(*vi, gh.graph()).s() << std::endl;
  }

  // bond parameters
  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = gh.edges(); ei != ei_end; ++ei) {
    std::cout << "bond " << *ei << ": type = " << gh.edge_type(*ei)
              << ", " << m.bond(*ei, gh.graph()) << std::endl;
  }
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::Parameters params;
  std::cin >> params;

  // get graph
  alps::graph_helper<> gh(params);

  // get model
  alps::model_helper<> mh(params);

  // construct from model library
  looper::model_parameter<> m0(params, gh.graph(), gh.inhomogeneous_sites(),
                               gh.inhomogeneous_bonds(), mh);
  output(gh, m0);

  // construct from parameters
  looper::model_parameter<> m1(gh.graph(), alps::half_integer<int>(1.5),
                               -2, -1);
  output(gh, m1);

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
