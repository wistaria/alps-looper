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

#include <looper/model.h>
#include <looper/util.h>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>

std::ostream& operator<<(std::ostream& os, const looper::bond_parameter& p)
{
  os << "C = " << p.c() << ", Jxy = " << p.jxy() << ", Jz = " << p.jz();
  return os;
}

std::ostream& operator<<(std::ostream& os, const looper::bond_matrix& m)
{
  boost::numeric::ublas::matrix<double> mat;
  looper::flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

int main()
{
  while (true) {
    double s0_in, s1_in, c, jxy, jz;
    std::cin >> s0_in >> s1_in >> c >> jxy >> jz;
    if (!std::cin) break;

    alps::half_integer<int> s0(s0_in);
    alps::half_integer<int> s1(s1_in);
    looper::bond_matrix xxz(s0, s1, looper::bond_parameter(c, jxy, jz));

    std::cout << "input parameters: S0 = " << s0 << ", S1 = " << s1
              << ", C = " << c << ", Jxy = " << jxy << ", Jz = " << jz
              << std::endl << xxz << std::endl;

    looper::bond_parameter p;
    bool success = looper::fit2bond(xxz.matrix(), p);

    assert(success);
    std::cout << "fitting result: " << p << std::endl;
    assert(p == looper::bond_parameter(c, jxy, jz));
  }
}
