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

#ifndef LOOPER_UTIL_H
#define LOOPER_UTIL_H

#include <boost/call_traits.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_float.hpp>
#include <boost/type_traits/is_integral.hpp>
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

template<typename T, typename U, typename R, typename A>
void flatten_matrix(const boost::multi_array<T, 2>& m_in,
                    boost::numeric::ublas::matrix<U, R, A>& m_out)
{
  assert(m_in.shape()[0] == m_in.shape()[1]);

  int dim = m_in.shape()[0];

  m_out.resize(dim, dim);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      m_out(i, j) = m_in[i][j];
}

template<typename T, typename U, typename R, typename A>
void flatten_matrix(const boost::multi_array<T, 4>& m_in,
                    boost::numeric::ublas::matrix<U, R, A>& m_out)
{
  assert(m_in.shape()[0] == m_in.shape()[2]);
  assert(m_in.shape()[1] == m_in.shape()[3]);

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
// function equal
//

#ifndef BOOST_NO_SFINAE

template<class T, class U>
bool equal(const T& x, const U& y,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0,
  typename boost::disable_if<boost::is_arithmetic<U> >::type* = 0)
{ return x == y; }

template<class T, class U>
bool equal(T x, U y,
  typename boost::enable_if<boost::is_integral<T> >::type* = 0,
  typename boost::enable_if<boost::is_integral<U> >::type* = 0)
{ return x == y; }

template<class T, class U>
bool equal(T x, U y, double tol,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{
  return (std::abs(x-y) < tol * std::abs(x)) ||
    (std::abs(x) < tol && std::abs(y) < tol);
}

template<class T, class U>
bool equal(T x, U y,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(x, y, 1.0e-12); }

template<class T, class U>
bool equal(const std::complex<T>& x, const std::complex<U>& y, double tol,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(x.real(), y.real(), tol) && equal(x.imag(), y.imag(), tol); }

template<class T, class U>
bool equal(const std::complex<T>& x, const std::complex<U>& y,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(x.real(), y.real()) && equal(x.imag(), y.imag()); }

template<class T, class U>
bool equal(const std::complex<T>& x, U y, double tol,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(x.real(), y, tol) && equal(x.imag(), 0., tol); }

template<class T, class U>
bool equal(const std::complex<T>& x, U y,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(x.real(), y) && equal(x.imag(), 0.); }

template<class T, class U>
bool equal(T x, const std::complex<U>& y, double tol,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(y, x, tol); }

template<class T, class U>
bool equal(T x, const std::complex<U>& y,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(y, x); }

#else

template<class T, class U>
bool equal(const T& x, const U& y)
{ return x == y; }

bool equal(double x, double y, double tol,
{
  return (std::abs(x-y) < tol * std::abs(x)) ||
    (std::abs(x) < tol && std::abs(y) < tol);
}

template<class T, class U>
bool equal(double x, double y,
  typename boost::enable_if<boost::is_float<T> >::type* = 0,
  typename boost::enable_if<boost::is_float<U> >::type* = 0)
{ return equal(x, y, 1.0e-12); }

template<class T, class U>
bool equal(const std::complex<T>& x, const std::complex<U>& y, double tol)
{ return equal(x.real(), y.real(), tol) && equal(x.imag(), y.imag(), tol); }

template<class T, class U>
bool equal(const std::complex<T>& x, const std::complex<U>& y)
{ return equal(x.real(), y.real()) && equal(x.imag(), y.imag()); }

template<class T>
bool equal(const std::complex<T>& x, double y, double tol)
{ return equal(x.real(), y, tol) && equal(x.imag(), 0., tol); }

template<class T>
bool equal(const std::complex<T>& x, double y)
{ return equal(x.real(), y) && equal(x.imag(), 0.); }

template<class T>
bool equal(double x, const std::complex<T>& y, double tol)
{ return equal(y, x, tol); }

template<class T>
bool equal(double x, const std::complex<T>& y)
{ return equal(y, x); }

#endif


//
// function is_zero
//

#ifndef BOOST_NO_SFINAE

template<class T>
bool is_zero(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return x == 0; }

template<class T>
bool is_zero(T x,
  typename boost::enable_if<boost::is_integral<T> >::type* = 0)
{ return x == 0; }

template<class T>
bool is_zero(T x, double tol,
  typename boost::enable_if<boost::is_float<T> >::type* = 0)
{ return std::abs(x) < tol; }

template<class T>
bool is_zero(T x,
  typename boost::enable_if<boost::is_float<T> >::type* = 0)
{ return is_zero(x, 1.0e-12); }

template<class T>
bool is_zero(const std::complex<T>& x, double tol)
{ return is_zero(x.real(), tol) && is_zero(x.imag(), tol); }

template<class T>
bool is_zero(const std::complex<T>& x)
{ return is_zero(x.real()) && is_zero(x.imag()); }

#else

template<class T>
bool is_zero(const T& x) { return x == 0; }

bool is_zero(double x, double tol) { return std::abs(x) < tol; }

bool is_zero(double x) { return is_zero(x, 1.0e-12); }

template<class T>
bool is_zero(const std::complex<T>& x, double tol)
{ return is_zero(x.real(), tol) && is_zero(x.imag(), tol); }

template<class T>
bool is_zero(const std::complex<T>& x)
{ return is_zero(x.real()) && is_zero(x.imag()); }

#endif


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
    if (!is_zero(x.imag()))
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

template<class T>
inline bool is_negative(const T& x,
  typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{ return x < 0; }
template<class T>
inline bool is_negative(T x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return x < 0; }
inline bool is_negative(unsigned char) { return false; }
inline bool is_negative(unsigned short) { return false; }
inline bool is_negative(unsigned int) { return false; }
inline bool is_negative(unsigned long) { return false; }

} // end namespace looper

#endif // LOOPER_UTIL_H
