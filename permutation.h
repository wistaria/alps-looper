// $Id: permutation.h 398 2003-10-09 10:33:05Z wistaria $

// looper-3 : C++ library for loop algorithm
//
// Copyright (C) 2001,2002  Synge Todo <wistaria@comp-phys.org>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

#ifndef LOOPER_ALGORITHM_PERMUTATION_H
#define LOOPER_ALGORITHM_PERMUTATION_H

// #include <looper/config.h>
#include <iterator>
#include <algorithm> // for std::iter_swap, std::random_shuffle

namespace looper {

template<class RandomAccessIter,
         class RandomNumberGenerator>
void random_shuffle(RandomAccessIter first, 
		    RandomAccessIter last,
		    RandomNumberGenerator& rng)
{
  typedef typename std::iterator_traits<RandomAccessIter>::difference_type
    diff_type;
  diff_type n = last - first;
  if (n < 2) return;

  for (std::size_t i = 1; i < n; ++i)
    std::iter_swap(first + i, first + boost::uniform_int<>(0, i)(rng));
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
  random_shuffle(perm_first, perm_first + c, rng);
  random_shuffle(perm_first + c, perm_last, rng);

  // reorder permutation according to values in guide0
  looper::guided_sort_binary(perm_first, perm_last,
			     guide0_first, guide0_last);
}

} // end namespace looper

#endif // LOOPER_PERMUTATION_H
