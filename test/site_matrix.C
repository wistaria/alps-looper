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

std::ostream& operator<<(std::ostream& os, const looper::site_matrix& m)
{
  boost::numeric::ublas::matrix<double> mat;
  looper::flatten_matrix(m.matrix(), mat);
  os << mat;
  return os;
}

std::ostream& operator<<(std::ostream& os, const looper::site_parameter& p)
{
  os << "C = " << p.c() << ", Hx = " << p.hx() << ", Hz = " << p.hz();
  return os;
}

int main()
{
  while (true) {
    double s_in, c, hx, hz;
    std::cin >> s_in >> c >> hx >> hz;
    if (!std::cin) break;

    looper::site_parameter s(s_in, c, hx, hz);
    looper::site_matrix site(s);

    std::cout << "input parameters: S = " << s.s()
              << ", C = " << c << ", Hx = " << hx << ", Hz = " << hz
              << std::endl << site << std::endl;

    looper::site_parameter p;
    bool success = looper::fit2site(site.matrix(), p);

    if (!success) std::cerr << "Error: fitting failed\n";

    std::cout << "fitting result: " << p << std::endl;
  }
}
