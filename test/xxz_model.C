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

#include <looper/xxz.h>
#include <iostream>

template<class G>
void output(const G& graph, const looper::xxz_model& m)
{
  typedef G graph_type;
  typedef typename boost::graph_traits<graph_type>::vertex_iterator
    vertex_iterator;
  typedef typename boost::graph_traits<graph_type>::edge_iterator
    edge_iterator;

  // site parameters
  std::cout << "number of spin types = " << m.num_spin_types() << std::endl;
  if (m.is_uniform_spin()) {
    std::cout << "S = " << m.uniform_spin() << std::endl;
  }
  typename alps::property_map<alps::site_type_t, graph_type, int>::const_type
    site_type(alps::get_or_default(alps::site_type_t(), graph, 0));
  vertex_iterator vi_end = boost::vertices(graph).second;
  for (vertex_iterator vi = boost::vertices(graph).first; vi != vi_end;
       ++vi) {
    int t = site_type[*vi];
    std::cout << "site " << *vi << ": type = " << t << ", S = " << m.spin(t)
              << std::endl;
  }

  // bond parameters
  std::cout << "number of bond types = " << m.num_bond_types() << std::endl;
  if (m.is_uniform_bond()) {
    std::cout << "C = " << m.uniform_bond().c()
              << ", Jxy = " << m.uniform_bond().jxy()
              << ", Jz = " << m.uniform_bond().jz()
              << std::endl;
  }
  typename alps::property_map<alps::bond_type_t, graph_type, int>::const_type
    bond_type(alps::get_or_default(alps::bond_type_t(), graph, 0));
  edge_iterator ei_end = boost::edges(graph).second;
  for (edge_iterator ei = boost::edges(graph).first; ei != ei_end; ++ei) {
    int t = bond_type[*ei];
    std::cout << "bond " << *ei << ": type = " << t
              << ", C = " << m.bond(t).c()
              << ", Jxy = " << m.bond(t).jxy()
              << ", Jz = " << m.bond(t).jz()
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
  typedef alps::graph_factory<>::graph_type graph_type;
  alps::graph_factory<graph_type> gf(params);
  const graph_type& graph = gf.graph();

  // get model
  alps::ModelLibrary models(params);

  // construct from model library
  looper::xxz_model m0(params, graph, models);
  output(graph, m0);

  // construct from parameters
  looper::xxz_model m1(-2, -1, alps::half_integer<short>(1.5), graph);
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
