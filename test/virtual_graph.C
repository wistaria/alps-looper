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

#include <looper/graph.h>
#include <alps/model.h>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef looper::graph_type rgraph_type;
  typedef boost::graph_traits<rgraph_type>::vertex_iterator rvertex_iterator;
  typedef boost::graph_traits<rgraph_type>::edge_iterator redge_iterator;

  typedef looper::virtual_graph<rgraph_type> vgraph_type;
  typedef boost::graph_traits<vgraph_type>::vertex_iterator vvertex_iterator;
  typedef boost::graph_traits<vgraph_type>::edge_iterator vedge_iterator;

  // real graph
  rgraph_type rg;
  looper::hypercubic_graph_generator<> gen(2, 2);
  looper::generate_graph(rg, gen);
  put(looper::vertex_type_t(), rg, *(vertices(rg).first), 1);
  set_parity(rg, looper::parity_t());
  std::cout << rg;
  rvertex_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi) {
    std::cout << looper::gauge(rg, *rvi) << ' ';
  }
  std::cout << std::endl;

  // virtual graph
  vgraph_type vg;
  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  vg.initialize(rg, spins);
  set_parity(vg);

  std::cout << "number of real vertices = "
            << num_vertices(rg) << std::endl;
  std::cout << "number of real edges = "
            << num_edges(rg) << std::endl;
  std::cout << "number of virtual vertices = "
            << num_vertices(vg) << std::endl;
  std::cout << "number of virtual edges = "
            << num_edges(vg) << std::endl;

  std::cout << vg;
  vvertex_iterator vvi, vvi_end;
  for (boost::tie(vvi, vvi_end) = vertices(vg); vvi != vvi_end; ++vvi) {
    std::cout << looper::gauge(vg, *vvi) << ' ';
  }
  std::cout << std::endl;
  vg.print_mapping(std::cout, rg);

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
