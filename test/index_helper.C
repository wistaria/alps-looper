/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: index_helper.C 554 2003-11-12 02:36:24Z wistaria $

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
  boost::mt19937 rng;
  boost::uniform_int<> uniform_int(1, 10);

  std::vector<int> a(100);
  check(a);
  std::cout << "std::vector(" << a.size() << ") OK\n";

  std::deque<int> b(uniform_int(rng));
  check(b);
  for (;;) {
    int n = uniform_int(rng);
    if (rng() < 0.5) {
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
