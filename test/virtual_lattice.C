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
    spin_wrapper(const value_type& v) : s(v) {}
    const value_type& s;
  };
  vector_spin_wrapper(const std::vector<value_type>& v) : vec_(v) {}
  template<class G>
  spin_wrapper site(const typename alps::graph_traits<G>::site_descriptor& v,
    const G& g) const
  { return spin_wrapper(vec_[boost::get(looper::site_type_t(), g, v)]); }
  const std::vector<value_type>& vec_;
};

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef looper::graph_type graph_type;
  typedef alps::graph_traits<graph_type>::site_iterator site_iterator;
  typedef alps::graph_traits<graph_type>::bond_iterator bond_iterator;

  // real graph
  graph_type rg;
  looper::hypercubic_graph_generator<> gen(2, 2);
  looper::generate_graph(rg, gen);
  put(looper::site_type_t(), rg, *(sites(rg).first), 1);
  set_parity(rg, looper::parity_t());
  std::cout << rg;
  site_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = sites(rg); rvi != rvi_end; ++rvi)
    std::cout << get(looper::parity_t(), rg, *rvi) << ' ';
  std::cout << std::endl;

  // virtual lattice
  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  looper::virtual_lattice<graph_type> vl(rg, vector_spin_wrapper<int>(spins));

  std::cout << "number of real sites = "
            << num_sites(rg) << std::endl;
  std::cout << "number of real bonds = "
            << num_bonds(rg) << std::endl;
  std::cout << "number of virtual sites = "
            << num_sites(vl) << std::endl;
  std::cout << "number of virtual bonds = "
            << num_bonds(vl) << std::endl;
  std::cout << vl;
  site_iterator vvi, vvi_end;
  for (boost::tie(vvi, vvi_end) = sites(vl); vvi != vvi_end; ++vvi)
    std::cout << gauge(vl, *vvi) << ' ';
  std::cout << std::endl;
  vl.mapping().output(std::cout, rg, vl.graph());

  vl.generate(rg, vector_spin_wrapper<int>(spins), true);

  std::cout << "number of real sites = "
            << num_sites(rg) << std::endl;
  std::cout << "number of real bonds = "
            << num_bonds(rg) << std::endl;
  std::cout << "number of virtual sites = "
            << num_sites(vl) << std::endl;
  std::cout << "number of virtual bonds = "
            << num_bonds(vl) << std::endl;
  std::cout << vl;
  for (boost::tie(vvi, vvi_end) = sites(vl); vvi != vvi_end; ++vvi)
    std::cout << gauge(vl, *vvi) << ' ';
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
