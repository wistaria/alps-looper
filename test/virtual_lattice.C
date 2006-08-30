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
  { return spin_wrapper(vec_[get(looper::site_type_t(), g, v)]); }
  const std::vector<value_type>& vec_;
};

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef looper::default_graph_type real_graph_type;
  typedef alps::graph_traits<real_graph_type>::site_iterator
    real_site_iterator;

  typedef looper::virtual_lattice<real_graph_type> virtual_lattice_type;
  typedef virtual_lattice_type::virtual_graph_type virtual_graph_type;
  typedef alps::graph_traits<virtual_graph_type>::site_iterator
    virtual_site_iterator;
  typedef alps::graph_traits<virtual_graph_type>::bond_iterator
    virtual_bond_iterator;

  // real graph
  alps::Parameters param;
  param["LATTICE"] = "square lattice";
  param["L"] = 2;
  alps::graph_helper<real_graph_type> gh(param);
  real_graph_type& g = gh.graph();
  put(looper::site_type_t(), g, *(sites(g).first), 1);
  std::cout << g;
  real_site_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = sites(g); rvi != rvi_end; ++rvi)
    std::cout << looper::gauge(gh, *rvi) << ' ';
  std::cout << std::endl;

  // virtual lattice
  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  virtual_lattice_type vl(gh, vector_spin_wrapper<int>(spins));

  std::cout << "number of real sites = "
            << num_sites(vl.rgraph()) << std::endl;
  std::cout << "number of real bonds = "
            << num_bonds(vl.rgraph()) << std::endl;
  std::cout << "number of virtual sites = "
            << num_vsites(vl) << std::endl;
  std::cout << "number of virtual bonds = "
            << num_vbonds(vl) << std::endl;
  std::cout << "maximum number of virutal sites = "
            << max_virtual_sites(vl) << std::endl;
  std::cout << vl;
  virtual_site_iterator vvi, vvi_end;
  for (boost::tie(vvi, vvi_end) = vsites(vl); vvi != vvi_end; ++vvi)
    std::cout << gauge(vl, *vvi) << ' ';
  std::cout << std::endl;
  vl.print_mapping(std::cout);
  for (boost::tie(vvi, vvi_end) = vsites(vl); vvi != vvi_end; ++vvi)
    std::cout << get(looper::site_index_t(), vl.rgraph(), rsite(vl, *vvi))
              << ' ';
  std::cout << std::endl;
  virtual_bond_iterator vei, vei_end;
  for (boost::tie(vei, vei_end) = vbonds(vl); vei != vei_end; ++vei)
    std::cout << get(looper::bond_index_t(), vl.rgraph(), rbond(vl, *vei))
              << ' ';
  std::cout << std::endl;

  vl.reinitialize(vector_spin_wrapper<int>(spins), true);

  std::cout << "number of real sites = "
            << num_sites(vl.rgraph()) << std::endl;
  std::cout << "number of real bonds = "
            << num_bonds(vl.rgraph()) << std::endl;
  std::cout << "number of virtual sites = "
            << num_vsites(vl) << std::endl;
  std::cout << "number of virtual bonds = "
            << num_vbonds(vl) << std::endl;
  std::cout << "maximum number of virutal sites = "
            << max_virtual_sites(vl) << std::endl;
  std::cout << vl;
  for (boost::tie(vvi, vvi_end) = vsites(vl); vvi != vvi_end; ++vvi)
    std::cout << gauge(vl, *vvi) << ' ';
  std::cout << std::endl;
  vl.print_mapping(std::cout);
  for (boost::tie(vvi, vvi_end) = vsites(vl); vvi != vvi_end; ++vvi)
    std::cout << get(looper::site_index_t(), vl.rgraph(), rsite(vl, *vvi))
              << ' ';
  std::cout << std::endl;
  for (boost::tie(vei, vei_end) = vbonds(vl); vei != vei_end; ++vei)
    if (is_real_bond(vl, *vei))
      std::cout << get(looper::bond_index_t(), vl.rgraph(), rbond(vl, *vei))
                << ' ';
  std::cout << std::endl;

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
