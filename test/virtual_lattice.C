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
  spin_wrapper site(
    const typename alps::graph_traits<G>::site_descriptor& v,
    const G& g) const
  {
    return spin_wrapper(vec_[get(looper::site_type_t(), g, v)]);
  }
  const std::vector<value_type>& vec_;
};

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef alps::coordinate_graph_type real_graph_type;
  typedef looper::lattice_helper<real_graph_type> lattice_type;
  typedef looper::real_site_iterator<lattice_type>::type real_site_iterator;
  typedef looper::real_bond_iterator<lattice_type>::type real_bond_iterator;
  typedef looper::virtual_site_iterator<lattice_type>::type
    virtual_site_iterator;
  typedef looper::virtual_bond_iterator<lattice_type>::type
    virtual_bond_iterator;

  alps::Parameters param;
  param["LATTICE"] = "square lattice";
  param["L"] = 2;
  lattice_type lat(param);
  put(looper::site_type_t(), const_cast<real_graph_type&>(lat.rg()), *(sites(lat.rg()).first), 1);
  std::cout << lat.rg();
  real_site_iterator rvi, rvi_end;
  for (boost::tie(rvi, rvi_end) = sites(lat.rg()); rvi != rvi_end; ++rvi)
    std::cout << looper::gauge(lat.rg(), *rvi) << ' ';
  std::cout << std::endl;

  std::vector<alps::half_integer<int> > spins(2);
  spins[0] = 1; spins[1] = 3./2;
  lat.generate_virtual_graph(vector_spin_wrapper<int>(spins));

  std::cout << "number of real sites = "
            << num_sites(lat.rg()) << std::endl;
  std::cout << "number of real bonds = "
            << num_bonds(lat.rg()) << std::endl;
  std::cout << "number of virtual sites = "
            << num_sites(lat.vg()) << std::endl;
  std::cout << "number of virtual bonds = "
            << num_bonds(lat.vg()) << std::endl;
  std::cout << "maximum number of virutal sites = "
            << max_virtual_sites(lat) << std::endl;
  std::cout << lat.vg();
  real_site_iterator rsi, rsi_end;
  for (boost::tie(rsi, rsi_end) = sites(lat.rg()); rsi != rsi_end; ++rsi) {
    virtual_site_iterator vsi, vsi_end;
    for (boost::tie(vsi, vsi_end) = sites(lat, *rsi); vsi != vsi_end;
         ++vsi)
      std::cout << get(looper::gauge_t(), lat.vg(), *vsi) << ' ';
  }
  std::cout << std::endl;
  lat.mp().output(std::cout, lat.rg(), lat.vg());
  for (boost::tie(rsi, rsi_end) = sites(lat.rg()); rsi != rsi_end; ++rsi) {
    virtual_site_iterator vsi, vsi_end;
    for (boost::tie(vsi, vsi_end) = sites(lat, *rsi); vsi != vsi_end;
         ++vsi)
      std::cout << get(looper::site_index_t(), lat.rg(),
                       get(looper::real_site_t(), lat.vg(), *vsi))
                << ' ';
  }
  std::cout << std::endl;
  real_bond_iterator rbi, rbi_end;
  for (boost::tie(rbi, rbi_end) = bonds(lat.rg()); rbi != rbi_end; ++rbi) {
    virtual_bond_iterator vbi, vbi_end;
    for (boost::tie(vbi, vbi_end) = bonds(lat, *rbi); vbi != vbi_end;
         ++vbi)
      std::cout << get(looper::bond_index_t(), lat.rg(),
                       get(looper::real_bond_t(), lat.vg(), *vbi))
                << ' ';
  }
  std::cout << std::endl;

  lat.generate_virtual_graph(vector_spin_wrapper<int>(spins), true);

  std::cout << "number of real sites = "
            << num_sites(lat.rg()) << std::endl;
  std::cout << "number of real bonds = "
            << num_bonds(lat.rg()) << std::endl;
  std::cout << "number of virtual sites = "
            << num_sites(lat.vg()) << std::endl;
  std::cout << "number of virtual bonds = "
            << num_bonds(lat.vg()) << std::endl;
  std::cout << "maximum number of virutal sites = "
            << max_virtual_sites(lat) << std::endl;
  std::cout << lat.vg();
  for (boost::tie(rsi, rsi_end) = sites(lat.rg()); rsi != rsi_end; ++rsi) {
    virtual_site_iterator vsi, vsi_end;
    for (boost::tie(vsi, vsi_end) = sites(lat, *rsi); vsi != vsi_end;
         ++vsi)
      std::cout << get(looper::gauge_t(), lat.vg(), *vsi) << ' ';
  }
  std::cout << std::endl;
  lat.mp().output(std::cout, lat.rg(), lat.vg());
  for (boost::tie(rsi, rsi_end) = sites(lat.rg()); rsi != rsi_end; ++rsi) {
    virtual_site_iterator vsi, vsi_end;
    for (boost::tie(vsi, vsi_end) = sites(lat, *rsi); vsi != vsi_end;
         ++vsi)
      std::cout << get(looper::site_index_t(), lat.rg(),
                       get(looper::real_site_t(), lat.vg(), *vsi))
                << ' ';
  }
  std::cout << std::endl;
  for (boost::tie(rbi, rbi_end) = bonds(lat.rg()); rbi != rbi_end; ++rbi) {
    virtual_bond_iterator vbi, vbi_end;
    for (boost::tie(vbi, vbi_end) = bonds(lat, *rbi); vbi != vbi_end;
         ++vbi)
      std::cout << get(looper::bond_index_t(), lat.rg(),
                       get(looper::real_bond_t(), lat.vg(), *vbi))
                << ' ';
  }
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
