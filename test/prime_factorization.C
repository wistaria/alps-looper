/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2010 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/expand.h>
#include <looper/prime_factorization.h>
#include <alps/utility/vectorio.hpp>
#include <numeric>
#include <cstdlib>
#include <iostream>

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  int n;
  while (true) {
    std::cin >> n;
    if (!std::cin) break;

    std::vector<int> factors = looper::prime_factorization(n);
    int r = std::accumulate(factors.begin(), factors.end(), 1, std::multiplies<int>());

    if (r == n)
      std::cout << n << " can be factored: " << n << " = ";
    else
      std::cout << n << " can not be factored: " << n << " < ";
    std::cout << alps::write_vector(factors, ":") << " = " << r << std::endl;
    if (r < n || r >= 2 * n) {
      std::cerr << "test fails\n";
      boost::throw_exception(std::logic_error("test fails"));
    }
  }

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
