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

/* $Id: permutation.h 693 2004-03-16 15:48:04Z wistaria $ */

#ifndef LOOPER_PERMUTATION_H
#define LOOPER_PERMUTATION_H

#include <iostream>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <algorithm> // for std::iter_swap, std::random_shuffle
#include <iterator>

namespace looper {

template<class RandomAccessIter,
         class RandomNumberGenerator>
void random_shuffle(RandomAccessIter first, 
                    RandomAccessIter last,
                    RandomNumberGenerator& rng)
{
  typedef typename std::iterator_traits<RandomAccessIter>::difference_type
    diff_type;

  for (diff_type i = 1; i < last - first; ++i) {
    std::iter_swap(first + i,
      first + boost::variate_generator<RandomNumberGenerator&,
        boost::uniform_int<> >(rng, boost::uniform_int<>(0, i))());
  }
}


// helper function : guided_sort_binary
//
// example
//   input  perm  : 0 1 2 3 4 5 6 7
//          guide : 1 0 0 1 0 1 1 0
//   output perm  : 7 1 2 4 3 5 6 0
//          guide : does not change

template<class RandomAccessIterator0, class RandomAccessIterator1>
void guided_sort_binary(RandomAccessIterator0 first_perm, 
                        RandomAccessIterator0 last_perm,
                        RandomAccessIterator1 first_guide,
                        RandomAccessIterator1 last_guide) {
  if (first_perm == last_perm) return;
  --last_perm;
  --last_guide;
  for (;;) {
    if (last_perm - first_perm <= 0) return;
    while (*first_guide == 0) {
      ++first_perm;
      ++first_guide;
      if (first_perm == last_perm) return;
    }
    while (*last_guide == 1) {
      --last_perm;
      --last_guide;
      if (first_perm == last_perm) return;
    }
    std::iter_swap(first_perm, last_perm);
    ++first_perm;
    ++first_guide;
    --last_perm;
    --last_guide;
  }
}


// function : restricted_random_shuffle
//
// This function generates a permutaion, which satisfies the
// conservation law.  This is used for generating a boundary graph for
// S >= 1 cases.  The model for RandomNumberGenerator is the same as
// for std::random_shuffle

template<class RandomAccessIter0, class RandomAccessIter1,
         class RandomNumberGenerator>
void restricted_random_shuffle(RandomAccessIter0 perm_first, 
                               RandomAccessIter0 perm_last,
                               RandomAccessIter1 guide0_first,
                               RandomAccessIter1 guide0_last,
                               RandomAccessIter1 guide1_first,
                               RandomAccessIter1 guide1_last,
                               RandomNumberGenerator& rng)
{
  typedef typename std::iterator_traits<RandomAccessIter1>::difference_type
    diff_type;

  // sort in two sectors (0 and 1) according to values in guide1
  guided_sort_binary(perm_first, perm_last, guide1_first, guide1_last);
  diff_type c = std::count(guide1_first, guide1_last, 0);

  // shuffle in each sector
  looper::random_shuffle(perm_first, perm_first + c, rng);
  looper::random_shuffle(perm_first + c, perm_last, rng);

  // reorder permutation according to values in guide0
  looper::guided_sort_binary(perm_first, perm_last,
                             guide0_first, guide0_last);
}

} // end namespace looper

#endif // LOOPER_PERMUTATION_H
