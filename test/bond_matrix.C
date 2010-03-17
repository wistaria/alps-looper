/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/flatten_matrix.h>
#include <looper/matrix.h>
#include <looper/model_parameter.h>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#ifdef HAVE_PARAPACK_13
# include <alps/math.hpp>
#else
# include <alps/numeric/is_equal.hpp>
#endif
#ifdef HAVE_PARAPACK_13
  using alps::is_equal;
#else
  using alps::numeric::is_equal;
#endif

std::ostream& operator<<(std::ostream& os, const looper::bond_matrix_xxz& m) {
  boost::numeric::ublas::matrix<double> mat;
  looper::flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

std::ostream& operator<<(std::ostream& os, const looper::bond_matrix_xyz& m) {
  boost::numeric::ublas::matrix<double> mat;
  looper::flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

int main() {
  while (true) {
    double s0_in, s1_in, c, jx, jy, jz;
    std::cin >> s0_in >> s1_in >> c >> jx >> jy >> jz;
    if (!std::cin) break;

    alps::half_integer<int> s0(s0_in);
    alps::half_integer<int> s1(s1_in);
    if (is_equal(jx, jy)) {
      looper::bond_parameter_xxz pi(c, jx, jz);
      looper::bond_matrix_xxz m(s0, s1, pi);
      std::cout << "[bond_matrix_xxz]\n";
      std::cout << "input parameters: S0 = " << s0 << ", S1 = " << s1 << ", "
                << pi << std::endl << m << std::endl;

      looper::bond_parameter_xxz pr(m.matrix());
      if (pi != pr) std::cerr << "Error: fitting failed\n";

      std::cout << "fitting result: " << pr << std::endl;
    }

    {
      std::cout << "[bond_matrix_xyz]\n";
      looper::bond_parameter_xyz pi(c, jx, jy, jz);
      looper::bond_matrix_xyz m(s0, s1, pi);
      std::cout << "input parameters: S0 = " << s0 << ", S1 = " << s1 << ", "
                << pi << std::endl << m << std::endl;

      looper::bond_parameter_xyz pr(m.matrix());
      if (pi != pr) std::cerr << "Error: fitting failed\n";

      std::cout << "fitting result: " << pr << std::endl;
    }
  }
}
