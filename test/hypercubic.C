/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/lattice.h>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

int main() {
  typedef alps::coordinate_graph_type graph_type;
  typedef looper::graph_traits<graph_type>::site_iterator site_iterator;

  std::vector<int> ext;
  int dim;
  std::cin >> dim;
  ext.resize(dim);
  for (std::vector<int>::iterator itr = ext.begin(); itr != ext.end(); ++itr) {
    std::cin >> *itr;
  }

  std::cout << "dimension = " << dim << std::endl;
  for (int d = 0; d < dim; ++d) {
    std::cout << "extent[" << d << "] = " << ext[d] << std::endl;
  }

  looper::hypercubic_graph_generator<> gen(ext);
  graph_type graph;
  looper::generate_graph(graph, gen);

  std::cout << graph;
}
