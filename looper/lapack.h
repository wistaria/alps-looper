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

#ifndef LOOPER_LAPACK_H
#define LOOPER_LAPACK_H

#include <algorithm>
#include <complex>
#include <stdexcept>

#include <alps/bindings/gels.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/lapack/heev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace {

template<class T> struct vector_helper;

template<class T>
struct vector_helper<boost::numeric::ublas::vector<T> >
{
  typedef T                                value_type;
  typedef boost::numeric::ublas::vector<T> vector_type;
  static value_type * begin_ptr(vector_type& v) { return &v(0); }
  static int size(const vector_type& v) { return v.size(); }
};

template<class T> struct matrix_helper;

template<class T, class A>
struct matrix_helper<
  boost::numeric::ublas::matrix<T,boost::numeric::ublas::row_major,A> >
{
  typedef T value_type;
  typedef boost::numeric::ublas::matrix<T,
            boost::numeric::ublas::row_major,A> matrix_type;
  BOOST_STATIC_CONSTANT(bool, is_row_major = true);
  BOOST_STATIC_CONSTANT(bool, is_column_major = false);
  static value_type * begin_ptr(matrix_type& m) { return &m(0, 0); }
  static int size1(const matrix_type& m) { return m.size1(); }
  static int size2(const matrix_type& m) { return m.size2(); }
};

template<class T, class A>
struct matrix_helper<
  boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major,A> >
{
  typedef T value_type;
  typedef boost::numeric::ublas::matrix<T,
            boost::numeric::ublas::column_major,A> matrix_type;
  BOOST_STATIC_CONSTANT(bool, is_row_major = false);
  BOOST_STATIC_CONSTANT(bool, is_column_major = true);
  static value_type * begin_ptr(matrix_type& m) { return &m(0, 0); }
  static int size1(const matrix_type& m) { return m.size1(); }
  static int size2(const matrix_type& m) { return m.size2(); }
};

}

namespace looper {

template <class Matrix, class Vector>
inline void diagonalize(Matrix& a, Vector& w, bool need_eigenvectors = true)
{
  typedef Matrix matrix_type;
  typedef Vector vector_type;

  // size check
  assert(matrix_helper<matrix_type>::size1(a) ==
         matrix_helper<matrix_type>::size2(a));
  assert(matrix_helper<matrix_type>::size1(a) ==
         vector_helper<vector_type>::size(w));

  char jobz;
  if (need_eigenvectors) {
    jobz = 'V';
  } else {
    jobz = 'N';
  }
  char uplo = 'L';
  int info;

  // call dispatcher
  info = boost::numeric::bindings::lapack::syev(jobz, uplo, a, w, boost::numeric::bindings::lapack::optimal_workspace());
  if (info != 0) throw std::runtime_error("failed in syev");
}

template <class T, class R, class A, class Vector>
inline void diagonalize(boost::numeric::ublas::matrix<std::complex<T>,R,A>& a, Vector& w, bool need_eigenvectors = true)
{
  typedef boost::numeric::ublas::matrix<std::complex<T>,R,A> matrix_type;
  typedef Vector                                             vector_type;

  // size check
  assert(matrix_helper<matrix_type>::size1(a) ==
         matrix_helper<matrix_type>::size2(a));
  assert(matrix_helper<matrix_type>::size1(a) ==
         vector_helper<vector_type>::size(w));

  char jobz;
  if (need_eigenvectors) {
    jobz = 'V';
  } else {
    jobz = 'N';
  }
  char uplo = 'L';
  int info;

  // call dispatcher
  info = boost::numeric::bindings::lapack::heev(jobz, uplo, a, w);
  if (info != 0) throw std::runtime_error("failed in heev");
}

template <class Matrix, class Vector>
inline void solve_llsp(Matrix& a, Vector& b, Vector& x)
{
  int m = matrix_helper<Matrix>::size1(a);
  int n = matrix_helper<Matrix>::size2(a);

  // size check
  assert(vector_helper<Vector>::size(b) == m);
  assert(vector_helper<Vector>::size(x) == n);

  char trans;
  int info;
  if (matrix_helper<Matrix>::is_column_major) {
    trans = 'N';
    if (m >= n) {
      info = boost::numeric::bindings::lapack::gels(trans, m, n, 1,
               matrix_helper<Matrix>::begin_ptr(a), m,
               vector_helper<Vector>::begin_ptr(b), m);
      std::copy(b.begin(), b.begin() + n, x.begin());
    } else {
      std::copy(b.begin(), b.begin() + m, x.begin());
      info = boost::numeric::bindings::lapack::gels(trans, m, n, 1,
               matrix_helper<Matrix>::begin_ptr(a), m,
               vector_helper<Vector>::begin_ptr(x), n);
    }
  } else {
    trans = 'T';
    if (m >= n) {
      info = boost::numeric::bindings::lapack::gels(trans, n, m, 1,
               matrix_helper<Matrix>::begin_ptr(a), n,
               vector_helper<Vector>::begin_ptr(b), m);
      std::copy(b.begin(), b.begin() + n, x.begin());
    } else {
      std::copy(b.begin(), b.begin() + m, x.begin());
      info = boost::numeric::bindings::lapack::gels(trans, n, m, 1,
               matrix_helper<Matrix>::begin_ptr(a), n,
               vector_helper<Vector>::begin_ptr(x), n);
    }
  }
  if (info != 0) throw std::runtime_error("failed in gels");
}

} // end namespace looper

#endif /* ! LOOPER_LAPACK_H */
