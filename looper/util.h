/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/spirit/core.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>
#include <limits>
#include <stdexcept>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

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
int alternating_tensor(boost::tuple<int, int, int> const& x)
{ return alternating_tensor(x.get<0>(), x.get<1>(), x.get<2>()); }


//
// function power2, power3, power4
//

#ifndef BOOST_NO_SFINAE

template<typename T>
T power2(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(std::complex<T> const& t)
{ return power2(real(t)) + power2(imag(t)); }

template<typename T>
T power3(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power3(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power4(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(std::complex<T> const& t)
{ return power2(power2(t)); }

#else

template<typename T>
T power2(T const& t) { return t * t; }

template<typename T>
T power3(T const& t) { return t * t * t; }

template<typename T>
T power4(T const& t) { return power2(power2(t)); }

#endif

//
// function flatten_matrix
//

template<typename T, typename U, typename R, typename A>
void flatten_matrix(boost::multi_array<T, 2> const& m_in,
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
void flatten_matrix(boost::multi_array<T, 4> const& m_in,
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
// dip (divide_if_positive)
//

template<typename T, typename U>
T dip(T x, U y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T x, U const& y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T const& x, U y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T const& x, U const& y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

//
// function crop_0, crop_01
//

template<typename T>
T crop_0(T x, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_0(T const& x, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_01(T x, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

template<typename T>
T crop_01(T const& x, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

//
// class tempalte integer_range
//

template<class T>
class integer_range
{
public:
  typedef T value_type;

  integer_range() {}
  integer_range(value_type v) : mi_(v), ma_(v) {}
  integer_range(integer_range const& r) : mi_(r.mi_), ma_(r.ma_) {}
  integer_range(std::string const& str,
                value_type def_mi = std::numeric_limits<value_type>::min(),
                value_type def_ma = std::numeric_limits<value_type>::max())
    : mi_(def_mi), ma_(def_ma)
  {
    using namespace boost::spirit;
    bool success;
    if (std::numeric_limits<value_type>::is_signed) {
      success = parse(str.c_str(),
        int_p[assign_a(mi_)][assign_a(ma_)] |
        ('[' >> !int_p[assign_a(mi_)] >> ':' >> !int_p[assign_a(ma_)] >> ']'),
        space_p).full;
    } else {
      success = parse(str.c_str(),
        uint_p[assign_a(mi_)][assign_a(ma_)] |
        ('[' >> !uint_p[assign_a(mi_)] >> ':' >> !uint_p[assign_a(ma_)] >> ']'),
        space_p).full;
    }
    if (!success) boost::throw_exception(std::runtime_error("parse error"));
  }

  value_type min() const { return mi_; }
  value_type max() const { return ma_; }

private:
  value_type mi_, ma_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class T>
std::ostream& operator<<(std::ostream& os, integer_range<T> const& ir)
{
  os << '[' << ir.min() << ':' << ir.max() << ']';
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

inline unsigned int generate_seed(int seed = -1)
{
  if (seed <= 0) {
    seed = boost::posix_time::microsec_clock::local_time().time_of_day().
      total_microseconds() << 24;
    seed &= ((1<<30)|((1<<30)-1));
  }
  return seed;
}

} // end namespace looper

#endif // LOOPER_UTIL_H
