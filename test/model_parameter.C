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

#include <looper/model.h>
#include <iostream>

template<typename G, typename SITE_P, typename BOND_P>
void output(const G& graph, const looper::model_parameter<SITE_P, BOND_P>& m)
{
  typedef G graph_type;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::edge_iterator
    edge_iterator;

  // site parameters
  std::cout << "number of site types = " << m.num_site_types() << std::endl;
  if (m.is_uniform_site()) {
    std::cout << "S = " << m.uniform_site().s() << std::endl;
  }
  typename alps::property_map<alps::site_type_t, graph_type, int>::const_type
    site_type(alps::get_or_default(alps::site_type_t(), graph, 0));
  vertex_iterator vi_end = boost::vertices(graph).second;
  for (vertex_iterator vi = boost::vertices(graph).first; vi != vi_end;
       ++vi) {
    int t = site_type[*vi];
    std::cout << "site " << *vi << ": type = " << t << ", S = "
              << m.site(t).s() << std::endl;
  }

  // bond parameters
  std::cout << "number of bond types = " << m.num_bond_types() << std::endl;
  if (m.is_uniform_bond()) std::cout << m.uniform_bond() << std::endl;
  typename alps::property_map<alps::bond_type_t, graph_type, int>::const_type
    bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));
  edge_iterator ei_end = boost::edges(graph).second;
  for (edge_iterator ei = boost::edges(graph).first; ei != ei_end; ++ei) {
    int t = bond_type[*ei];
    std::cout << "bond " << *ei << ": type = " << t << ", " << m.bond(t)
              << std::endl;
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
  typedef alps::graph_helper<>::graph_type graph_type;
  alps::graph_helper<> gh(params);
  const graph_type& graph = gh.graph();

  // get model
  alps::ModelLibrary models(params);

  // construct from model library
  looper::model_parameter<> m0(params, graph, models, false);
  output(graph, m0);

  // construct from parameters
  looper::model_parameter<> m1(-2, -1, alps::half_integer<int>(1.5), graph, false);
  output(graph, m1);

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