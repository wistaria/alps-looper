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

/* $Id: index_helper.C 698 2004-03-17 09:23:58Z wistaria $ */

#include <looper/vector_helper.h>
#include <boost/random.hpp>
#include <iostream>

template<class C>
void check(const C& c)
{
  looper::index_helper<C> helper(c);
  for (int i = 0; i < c.size(); ++i) assert(i == helper.index(&c[i]));
}

int main()
{
  // random number generators
  boost::mt19937 base_rng;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(base_rng, boost::uniform_real<>());
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
    uniform_int(base_rng ,boost::uniform_int<>(1, 10));

  std::vector<int> a(100);
  check(a);
  std::cout << "std::vector(" << a.size() << ") OK\n";

  std::deque<int> b(uniform_int());
  check(b);
  for (;;) {
    int n = uniform_int();
    if (uniform_01() < 0.5) {
      for (int i = 0; i < n; ++i) b.push_front(0);
    } else {
      for (int i = 0; i < n; ++i) b.push_back(0);
    }
    check(b);
    std::cout << "std::deque(" << b.size() << ") OK\n";
    if (b.size() > 300) break;
  }
  assert(looper::index_helper<>::index(b, &b[b.size() / 2]) == b.size() / 2);
}
