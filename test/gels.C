/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <iostream>

int main()
{
  // random number generator
  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    rng(boost::mt19937(), boost::uniform_real<>(-.5, .5));
  for (int i = 0; i < 1000; ++i) rng();

  int m, n;
  std::cin >> m >> n;

  std::cout << "m: dimension of problem    : " << m << std::endl;
  std::cout << "n: number of basis vectors : " << n << std::endl;

  boost::numeric::ublas::matrix<double> a(m, n);
  boost::numeric::ublas::vector<double> b(m);
  boost::numeric::ublas::vector<double> x(n);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) a(i,j) = rng();
    b(i) = rng();
  }

  std::cout << "input:  A   = " << a << std::endl;
  std::cout << "input:  B   = " << b << std::endl;

  boost::numeric::ublas::matrix<double> atmp(a);
  looper::solve_llsp(atmp, b, x);

  std::cout << "output: X   = " << x << std::endl;
  std::cout << "check:  AxX = " << prod(a, x) << std::endl;
}
