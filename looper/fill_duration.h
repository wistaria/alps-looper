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

#ifndef LOOPER_FILL_DURATION_H
#define LOOPER_FILL_DURATION_H

#include <iostream>
#include <boost/random/exponential_distribution.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

// fill duration [0,tmax] uniformly with density r

template<class RNG, class C>
void fill_duration(RNG& uniform_01, C& array,
                   typename C::value_type r, typename C::value_type tmax)
{
  typedef typename C::value_type value_type;

  if (tmax < value_type(0.))
    boost::throw_exception(std::invalid_argument("invalid argument"));

  array.clear();
  if (r <= value_type(0.)) return;

  boost::exponential_distribution<value_type> exp_rng(r);
  value_type t = value_type(0.);
  while (true) {
    t += exp_rng(uniform_01);
    if (t >= tmax) break;
    array.push_back(t);
  }
}

// fill duration [0,1] uniformly with density r

template<class RNG, class C>
void fill_duration(RNG& uniform_01, C& array, typename C::value_type r)
{
  typedef typename C::value_type value_type;
  fill_duration(uniform_01, array, r, value_type(1.));
}

} // end namespace looper

#endif // LOOPER_FILL_DURATION_H
