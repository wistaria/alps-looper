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
#include <boost/numeric/ublas/io.hpp>
#include <boost/random.hpp>
#include <iostream>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace boost::numeric::ublas;
#endif

int main()
{
  const double tol = 1.0e-10;

  // matrix & vector
  typedef boost::numeric::ublas::matrix<double,
    boost::numeric::ublas::column_major> matrix_type;
  typedef boost::numeric::ublas::vector<double> vector_type;

  // random number generator
  boost::variate_generator<boost::mt19937, boost::uniform_real<> >
    rng(boost::mt19937(4357), boost::uniform_real<>(-.5, .5));

  while (true) {
    int m, n;
    std::cin >> m >> n;
    if (!std::cin) break;

    std::cout << "m: dimension of problem    : " << m << std::endl;
    std::cout << "n: number of basis vectors : " << n << std::endl;

    matrix_type a(m, n);
    vector_type b(m);
    vector_type x(n);
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) a(i,j) = rng();
      b(i) = rng();
    }

    std::cout << "input: A      = " << a << std::endl;
    std::cout << "input: B      = " << b << std::endl;

    double r = looper::solve_llsp(a, b, x, tol);
    if (r < tol) r = 0.0;

    vector_type bt = prod(a,x);
    std::cout << "output: X     = " << x << std::endl;
    std::cout << "check: A*X    = " << bt << std::endl
              << "       B      = " << b << std::endl;
    std::cout << "residual      = " << r << std::endl;
  }
}
