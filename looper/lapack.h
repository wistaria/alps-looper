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

#include <looper/config.h>

#include <algorithm>
#include <complex>
#include <stdexcept>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#ifdef HAVE_BINDINGS

#include <alps/bindings/gels.hpp>
#include <alps/bindings/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>

#else

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

void dsyev_ (const char& jobz, const char& uplo, const int& n,
             double a[], const int& lda, double w[], double work[], const int& lwork,
             int& info);


void dgels_ (const char& trans, const int& m, const int& n, const int& nrhs,
             double a[], const int& lda, double b[], const int& ldb,
             double work[], const int& lwork, int& info);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

namespace looper {

namespace lapack {

inline void syev(const char& jobz, const char& uplo, const int& n,
                 double a[], const int& lda, double w[],
                 int& info) {
  // check optimal size of lwork
  int lwork = -1;
  double* work  = new double[1];
  dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
  lwork = static_cast<int>(work[0]);
  delete [] work;

  // do real work
  work  = new double[lwork];
  dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
  delete [] work;
}

inline void gels(const char& trans, const int& m, const int& n,
                 const int& nrhs,
                 double a[], const int& lda,
                 double b[], const int& ldb,
                 int& info) {
  int lwork = std::min(m,n) + std::max(std::max(1,m), std::max(n,nrhs)) * 64;
  double* work  = new double[lwork];
  dgels_(trans, m, n, nrhs, a, lda, b, ldb, work, lwork,info);
  delete [] work;
}

} // namespace lapack

} // namespace looper

#endif

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
  // size check
  assert(matrix_helper<Matrix>::size1(a) == matrix_helper<Matrix>::size2(a));
  assert(matrix_helper<Matrix>::size1(a) == vector_helper<Vector>::size(w));

  char jobz;
  if (need_eigenvectors) {
    jobz = 'V';
  } else {
    jobz = 'N';
  }
  char uplo = 'L';
  int info;

  // call dispatcher
#ifdef HAVE_BINDINGS
  info = boost::numeric::bindings::lapack::syev(jobz, uplo, a, w);
#else
  lapack::syev(jobz, uplo,
               vector_helper<Vector>::size(w),
               matrix_helper<Matrix>::begin_ptr(a),
               vector_helper<Vector>::size(w),
               vector_helper<Vector>::begin_ptr(w),
               info);
#endif
  if (info != 0) throw std::runtime_error("failed in syev");
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
#ifdef HAVE_BINDINGS
      info = boost::numeric::bindings::lapack::gels(trans, m, n, 1,
               matrix_helper<Matrix>::begin_ptr(a), m,
               vector_helper<Vector>::begin_ptr(b), m);
#else
      lapack::gels(trans, m, n, 1,
                   matrix_helper<Matrix>::begin_ptr(a), m,
                   vector_helper<Vector>::begin_ptr(b), m,
                   info);
#endif
      std::copy(b.begin(), b.begin() + n, x.begin());
    } else {
      std::copy(b.begin(), b.begin() + m, x.begin());
#ifdef HAVE_BINDINGS
      info = boost::numeric::bindings::lapack::gels(trans, m, n, 1,
               matrix_helper<Matrix>::begin_ptr(a), m,
               vector_helper<Vector>::begin_ptr(x), n);
#else
      lapack::gels(trans, m, n, 1,
                   matrix_helper<Matrix>::begin_ptr(a), m,
                   vector_helper<Vector>::begin_ptr(x), n,
                   info);
#endif
    }
  } else {
    trans = 'T';
    if (m >= n) {
#ifdef HAVE_BINDINGS
      info = boost::numeric::bindings::lapack::gels(trans, n, m, 1,
               matrix_helper<Matrix>::begin_ptr(a), n,
               vector_helper<Vector>::begin_ptr(b), m);
#else
      lapack::gels(trans, n, m, 1,
                   matrix_helper<Matrix>::begin_ptr(a), n,
                   vector_helper<Vector>::begin_ptr(b), m,
                   info);
#endif
      std::copy(b.begin(), b.begin() + n, x.begin());
    } else {
      std::copy(b.begin(), b.begin() + m, x.begin());
#ifdef HAVE_BINDINGS
      info = boost::numeric::bindings::lapack::gels(trans, n, m, 1,
               matrix_helper<Matrix>::begin_ptr(a), n,
               vector_helper<Vector>::begin_ptr(x), n);
#else
      lapack::gels(trans, n, m, 1,
                   matrix_helper<Matrix>::begin_ptr(a), n,
                   vector_helper<Vector>::begin_ptr(x), n,
                   info);
#endif
    }
  }
  if (info != 0) throw std::runtime_error("failed in gels");
}

} // end namespace looper

#endif // LOOPER_LAPACK_H
