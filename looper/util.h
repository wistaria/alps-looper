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

#ifndef LOOPER_UTIL_H
#define LOOPER_UTIL_H

#include <boost/call_traits.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <complex>
#include <stdexcept>

namespace looper {

//
// function alternating_tensor
//

inline int alternating_tensor(int i, int j, int k)
{
  switch (i) {
  case 0 :
    switch (j) {
    case 1:
      return (k == 2) ? 1 : 0;
    case 2:
      return (k == 1) ? -1 : 0;
    default:
      return 0;
    }
  case 1 :
    switch (j) {
    case 2:
      return (k == 0) ? 1 : 0;
    case 0:
      return (k == 2) ? -1 : 0;
    default:
      return 0;
    }
  case 2 :
    switch (j) {
    case 0:
      return (k == 1) ? 1 : 0;
    case 1:
      return (k == 0) ? -1 : 0;
    default:
      return 0;
    }
  default :
    break;
  }
  return 0;
}

inline
int alternating_tensor(const boost::tuple<int, int, int>& x)
{
  return alternating_tensor(x.get<0>(), x.get<1>(), x.get<2>());
}


//
// function sqr
//

#ifndef BOOST_NO_SFINAE

template<typename T>
T sqr(T t, typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{
  return t * t;
}

template<typename T>
T sqr(const T& t, typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{
  return t * t;
}

#else

template<typename T>
T sqr(const T& t) { return t * t; }

#endif

//
// function flatten_matrix
//

template<typename T, typename U>
void flatten_matrix(const boost::multi_array<T, 4>& m_in,
                    boost::numeric::ublas::matrix<U>& m_out)
{
#ifndef NDEBUG
  assert(m_in.shape()[0] == m_in.shape()[2]);
  assert(m_in.shape()[1] == m_in.shape()[3]);
#endif

  int d0 = m_in.shape()[0];
  int d1 = m_in.shape()[1];
  int dim = d0 * d1;

  m_out.resize(dim, dim);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
        for (int j1 = 0; j1 < d1; ++j1)
          m_out(i0 * d1 + i1, j0 * d1 + j1) = m_in[i0][i1][j0][j1];
}


//
// function nearly_equal
//

inline bool nearly_equal(double x, double y, double tol = 1.0e-12)
{
  return (std::abs(x-y) < tol * std::abs(x)) ||
    (std::abs(x) < tol && std::abs(y) < tol);
}


//
// function nearly_zero
//

inline bool nearly_zero(double x, double tol = 1.0e-12)
{
  return std::abs(x) < tol;
}


//
// function numeric_cast
//

namespace detail {

template<typename U, typename T>
struct numeric_cast_helper {
  static U value(typename boost::call_traits<T>::param_type x)
  {
    return x;
  }
};

template<typename U, typename T>
struct numeric_cast_helper<U, std::complex<T> > {
  static U value(const std::complex<T>& x) {
    if (x.imag() != 0)
      boost::throw_exception(std::runtime_error("can not convert complex number into real one"));
    return x.real();
  }
};

} // end namespace detail


//
// function range_01
//

#ifndef BOOST_NO_SFINAE

template<typename T>
T range_01(T x, typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{
  typedef T value_type;
  using std::min; using std::max;
  return min(max(x, value_type(0)), value_type(1));
}

template<typename T>
T range_01(const T& x, typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{
  typedef T value_type;
  using std::min; using std::max;
  return min(max(x, value_type(0)), value_type(1));
}

#else

template<typename T>
T range_01(const T& x)
{
  typedef T value_type;
  using std::min; using std::max;
  return min(max(x, value_type(0)), value_type(1));
}

#endif

} // end namespace looper

#endif // LOOPER_UTIL_H
