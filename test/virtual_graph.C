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

#include <looper/lattice.h>
#include <alps/model.h>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace alps;
#endif

template<class I>
struct vector_spin_wrapper
{
  typedef alps::half_integer<I> value_type;
  struct spin_wrapper
  {
    spin_wrapper(const value_type& v) : val_(v) {}
    const value_type& s() const { return val_; }
    const value_type& val_;
  };
  vector_spin_wrapper(const std::vector<value_type>& v) : vec_(v) {}
  template<class G>
  spin_wrapper site(const typename boost::graph_traits<G>::vertex_descriptor& v,
    const G& g) const
  { return spin_wrapper(vec_[boost::get(looper::vertex_type_t(), g, v)]); }
  const std::vector<value_type>& vec_;
};

struct unity_weight
{
  template<class G>
  double operator()(typename boost::graph_traits<G>::vertex_descriptor,
    const G&) const
  { return 1; }
  template<class G>
  double operator()(typename boost::graph_traits<G>::edge_descriptor,
    const G&) const
  { return 1; }
};

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef looper::graph_type graph_type;
  typedef boost::graph_traits<graph_type>::vertex_iterator vertex_iterator;
  typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;

  // real graph
  graph_type rg;
  looper::hypercubic_graph_generator<> gen(2, 2);
  looper::generate_graph(rg, gen);
  put(looper::vertex_type_t(), rg, *(vertices(rg).first), 1);
  set_parity(rg, looper::parity_t());
  std::cout << rg;
  vertex_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = vertices(rg); rvi != rvi_end; ++rvi) {
    std::cout << get(looper::parity_t(), rg, *rvi) << ' ';
  }
  std::cout << std::endl;

  // virtual lattice
  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  // graph_type vg;
  // looper::virtual_mapping<graph_type> vm;
  // looper::generate_virtual_lattice(vg, vm, rg,
  //   vector_spin_wrapper<int>(spins), unity_weight());
  looper::virtual_lattice<graph_type> vl(rg, vector_spin_wrapper<int>(spins), unity_weight());
  // set_parity(vg, looper::parity_t());
  set_parity(vl, looper::parity_t());

  std::cout << "number of real vertices = "
            << num_vertices(rg) << std::endl;
  std::cout << "number of real edges = "
            << num_edges(rg) << std::endl;
  std::cout << "number of virtual vertices = "
            << num_vertices(vl) << std::endl;
  std::cout << "number of virtual edges = "
            << num_edges(vl) << std::endl;

  std::cout << vl;
  vertex_iterator vvi, vvi_end;
  for (boost::tie(vvi, vvi_end) = vertices(vl); vvi != vvi_end; ++vvi) {
    // std::cout << get(looper::parity_t(), vl.graph(), *vvi) << ' ';
    std::cout << gauge(vl, *vvi) << ' ';
  }
  std::cout << std::endl;
  vl.mapping().output(std::cout, rg, vl.graph());

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
