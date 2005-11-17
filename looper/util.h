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
#include <boost/spirit/core.hpp>
#include <boost/throw_exception.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_float.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>
#include <limits>
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
// dip (divide_if_positive)
//

template<typename T>
T dip(T x, T y, typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{ return (y > T(0)) ? (x / y) : T(0); }

template<typename T>
T dip(const T& x, const T& y,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return (y > T(0)) ? (x / y) : T(0); }


//
// function crop_0, crop_01
//

template<typename T>
T crop_0(T x, typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_0(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_01(T x, typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

template<typename T>
T crop_01(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

//
// class tempalte integer_range
//

template<class T>
class integer_range {
public:
  typedef T value_type;

  integer_range() : mi_(), ma_() {}
  integer_range(value_type v) : mi_(v), ma_(v) {}
  integer_range(const integer_range& r) : mi_(r.mi_), ma_(r.ma_) {}
  integer_range(const std::string& str,
                value_type def_mi = std::numeric_limits<value_type>::min(),
                value_type def_ma = std::numeric_limits<value_type>::max()) :
    mi_(), ma_()
  {
    using namespace boost::spirit;
    bool success;
    value_type mi = def_mi;
    value_type ma = def_ma;
    if (std::numeric_limits<value_type>::is_signed) {
      success = parse(str.c_str(),
        int_p[assign_a(mi)][assign_a(ma)] |
        ('[' >> !int_p[assign_a(mi)] >> ':' >> !int_p[assign_a(ma)] >> ']'),
        space_p).full;
    } else {
      success = parse(str.c_str(),
        uint_p[assign_a(mi)][assign_a(ma)] |
        ('[' >> !uint_p[assign_a(mi)] >> ':' >> !uint_p[assign_a(ma)] >> ']'),
        space_p).full;
    }
    if (!success)
      boost::throw_exception(std::runtime_error("parse error"));
    mi_ = mi;
    ma_ = ma;
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
std::ostream& operator<<(std::ostream& os, const integer_range<T>& ir)
{
  os << '[' << ir.min() << ':' << ir.max() << ']';
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_UTIL_H
