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

#include <looper/union_find.h>
#include <cmath>
#include <boost/random.hpp>
#include <iostream>
#include <valarray>

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
using namespace looper::union_find;
#endif

const int n = 100;

template<class Itr0, class Itr1>
int index(const Itr0& itr, const Itr1& base) {
  return &(*itr) - &(*base);
}

int main()
{
  typedef std::vector<looper::union_find::node<> > vector_type;

  // random number generator
  boost::variate_generator<boost::mt19937, boost::uniform_int<> >
    rng(boost::mt19937(4357), boost::uniform_int<>(0, n - 1));

  std::cout << "[[union find test]]\n";

  vector_type tree(n);

  std::cout << "\n[making tree]\n";

  for (int i = 0; i < n; i++) {
    int i0 = rng();
    int i1 = rng();
    std::cout << "connecting node " << i0 << " to node " << i1 << std::endl;
    unify(tree[i0], tree[i1]);
  }

  std::cout << "\n[results]\n";

  for (vector_type::iterator itr = tree.begin(); itr != tree.end(); ++itr) {
    if (itr->is_root()) {
      std::cout << "node " << index(itr, tree.begin())
                << " is root and tree size is "
                << itr->weight() << std::endl;
    } else {
      std::cout << "node " << index(itr, tree.begin());
      std::cout << "'s parent is " << index(itr->parent(), tree.begin());
      std::cout << " and its root is " << index(itr->root(), tree.begin())
                << std::endl;
    }
  }
}
