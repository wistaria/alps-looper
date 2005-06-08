/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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
  // random number generator
  boost::mt19937 eng(29833u);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
    rng(eng, boost::uniform_int<>(0, n-1));

  std::cout << "[[union find test]]\n";

  std::vector<looper::union_find::node_idx> nodes_idx(n);
  std::vector<looper::union_find::node<> > nodes_ptr(n);

  std::cout << "\n[making tree]\n";

  for (int i = 0; i < n; i++) {
    int i0 = rng();
    int i1 = rng();
    std::cout << "connecting node " << i0 << " to node " << i1 << std::endl;
    looper::union_find::unify(nodes_idx, i0, i1);
    looper::union_find::unify(nodes_ptr[i0], nodes_ptr[i1]);
  }

  std::cout << "\n[results (index based)]\n";

  for (std::vector<looper::union_find::node_idx>::iterator
	 itr = nodes_idx.begin(); itr != nodes_idx.end(); ++itr) {
    int g = index(itr, nodes_idx.begin());
    if (itr->is_root()) {
      std::cout << "node " << g
		<< " is root and tree size is "
                << itr->weight() << std::endl;
    } else {
      std::cout << "node " << g
		<< "'s parent is " << itr->parent
		<< " and its root is " << root_index(nodes_idx, g)
                << std::endl;
    }
  }

  std::cout << "\n[results (pointer based)]\n";

  for (std::vector<looper::union_find::node<> >::iterator itr = nodes_ptr.begin();
       itr != nodes_ptr.end(); ++itr) {
    if (itr->is_root()) {
      std::cout << "node " << index(itr, nodes_ptr.begin())
                << " is root and nodes_ptr size is "
                << itr->weight() << std::endl;
    } else {
      std::cout << "node " << index(itr, nodes_ptr.begin())
		<< "'s parent is " << index(itr->parent(), nodes_ptr.begin())
		<< " and its root is " << index(itr->root(), nodes_ptr.begin())
                << std::endl;
    }
  }
}
