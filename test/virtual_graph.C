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

/* $Id: virtual_graph.C 717 2004-03-23 09:16:54Z wistaria $ */

#include <looper/graph.h>
#include <looper/virtual_graph.h>
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
  looper::virtual_graph<graph_type> vg;
  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  looper::generate_virtual_graph(vg, rg, spins);
  alps::set_parity(vg);

  std::cout << "number of original real vertices = " << vg.num_real_vertices
            << std::endl;
  std::cout << "number of original real edges = " << vg.num_real_edges
            << std::endl;
  std::cout << vg.graph;
  for (vertex_iterator vi = boost::vertices(vg.graph).first;
       vi != boost::vertices(vg.graph).second; ++vi) {
    std::cout << looper::gauge(*vi, vg.graph) << ' ';
    //    std::cout << (boost::get(looper::parity_t(), vg.graph, *vi) == alps::parity::black ? 0 : 1) << ' ';
  }
  std::cout << std::endl;
  std::cout << vg.mapping;

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
