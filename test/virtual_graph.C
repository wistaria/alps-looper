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

  typedef looper::parity_graph_type graph_type;
  typedef graph_type::vertex_iterator vertex_iterator;
  typedef graph_type::edge_iterator edge_iterator;

  // real graph
  looper::hypercubic_graph_generator<> gen(2, 2);
  graph_type rg;
  looper::generate_graph(rg, gen);
  boost::put(looper::vertex_type_t(), rg, *(boost::vertices(rg).first), 1);
  alps::set_parity(rg);
  std::cout << rg;
  for (vertex_iterator vi = boost::vertices(rg).first;
       vi != boost::vertices(rg).second; ++vi) {
    std::cout << looper::gauge(*vi, rg) << ' ';
  }
  std::cout << std::endl;

  // virtual graph
  graph_type vg;
  looper::virtual_mapping<graph_type> vm;
  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  looper::generate_virtual_graph(rg, spins, vg, vm);
  alps::set_parity(vg);

  std::cout << "number of real vertices = "
            << boost::num_vertices(rg) << std::endl;
  std::cout << "number of real edges = "
            << boost::num_edges(rg) << std::endl;
  std::cout << "number of virtual vertices = "
            << boost::num_vertices(vg) << std::endl;
  std::cout << "number of virtual edges = "
            << boost::num_edges(vg) << std::endl;

  std::cout << vg;
  for (vertex_iterator vi = boost::vertices(vg).first;
       vi != boost::vertices(vg).second; ++vi) {
    std::cout << looper::gauge(*vi, vg) << ' ';
  }
  std::cout << std::endl;
  vm.output(std::cout, rg, vg);

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
