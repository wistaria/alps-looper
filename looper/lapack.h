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

#include <looper/fortran.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <algorithm>
#include <complex>
#include <stdexcept>

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

//
// Simple Driver Routines for Eigenvalue Problems for Symmetric/Hermitian
// Matrixes
//

// ?SYEV: Computes all eigenvalues and, optionally, eigenvectors of a
// real symmetric/complex Hermite matrix

void LOOPER_FCALL(ssyev,SSYEV)(const char& jobz, const char& uplo,
                               const int& n,
                               float a[], const int& lda,
                               float w[],
                               float work[], const int& lwork,
                               int& info);

void LOOPER_FCALL(dsyev,DSYEV)(const char& jobz, const char& uplo,
                               const int& n,
                               double a[], const int& lda,
                               double w[],
                               double work[], const int& lwork,
                               int& info);

void LOOPER_FCALL(cheev,CHEEV)(const char& jobz, const char& uplo,
                               const int& n,
                               std::complex<float> a[], const int& lda,
                               float w[],
                               std::complex<float> work[], const int& lwork,
                               float rwork[],
                               int& info);

void LOOPER_FCALL(zheev,ZHEEV)(const char& jobz, const char& uplo,
                               const int& n,
                               std::complex<double> a[], const int& lda,
                               double w[],
                               std::complex<double> work[], const int& lwork,
                               double rwork[],
                               int& info);

//
// Simple Driver Routines for Solving Least Squares Problems
//

// ?GELS: Computes the least squares solution to an over-determined system
// of linear equations, A X=B or A**H X=B, or the minimum norm solution of
// an under-determined system, where A is a general rectangular matrix of
// full rank, using a QR or LQ factorization of A.

void LOOPER_FCALL(sgels,SGELS)(const char& trans, const int& m, const int& n,
                               const int& nrhs,
                               float a[], const int& lda,
                               float b[], const int& ldb,
                               float work[], const int& lwork,
                               int& info);

void LOOPER_FCALL(dgels,DGELS)(const char& trans, const int& m, const int& n,
                               const int& nrhs,
                               double a[], const int& lda,
                               double b[], const int& ldb,
                               double work[], const int& lwork,
                               int& info);

void LOOPER_FCALL(cgels,CGELS)(const char& trans, const int& m, const int& n,
                               const int& nrhs,
                               std::complex<float> a[], const int& lda,
                               std::complex<float> b[], const int& ldb,
                               std::complex<float> work[], const int& lwork,
                               int& info);

void LOOPER_FCALL(zgels,ZGELS)(const char& trans, const int& m, const int& n,
                               const int& nrhs,
                               std::complex<double> a[], const int& lda,
                               std::complex<double> b[], const int& ldb,
                               std::complex<double> work[], const int& lwork,
                               int& info);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

namespace looper {

namespace lapack_dispatch {

inline void syev(const char& jobz, const char& uplo, const int& n,
                 float a[], const int& lda, float w[],
                 int& info) {
  // check optimal size of lwork
  int lwork = -1;
  float* work  = new float[1];
  LOOPER_FCALL(ssyev,SSYEV)(jobz, uplo, n, a, lda, w, work, lwork, info);
  lwork = static_cast<int>(work[0]);
  delete [] work;

  // do real work
  work  = new float[lwork];
  LOOPER_FCALL(ssyev,SSYEV)(jobz, uplo, n, a, lda, w, work, lwork, info);
  delete [] work;
}

inline void syev(const char& jobz, const char& uplo, const int& n,
                 double a[], const int& lda, double w[],
                 int& info) {
  // check optimal size of lwork
  int lwork = -1;
  double* work  = new double[1];
  LOOPER_FCALL(dsyev,DSYEV)(jobz, uplo, n, a, lda, w, work, lwork, info);
  lwork = static_cast<int>(work[0]);
  delete [] work;

  // do real work
  work  = new double[lwork];
  LOOPER_FCALL(dsyev,DSYEV)(jobz, uplo, n, a, lda, w, work, lwork, info);
  delete [] work;
}

inline void syev(const char& jobz, const char& uplo, const int& n,
                 std::complex<float> a[], const int& lda, float w[],
                 int& info) {
  // check optimal size of lwork
  int lwork = -1;
  std::complex<float>* work  = new std::complex<float>[1];
  float* rwork = new float[std::max(1, 3*n-2)];
  LOOPER_FCALL(cheev,CHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  lwork = static_cast<int>(std::real(work[0]));
  delete [] work;

  // do real work
  work  = new std::complex<float>[lwork];
  LOOPER_FCALL(cheev,CHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  delete [] rwork;
  delete [] work;
}

inline void syev(const char& jobz, const char& uplo, const int& n,
                 std::complex<double> a[], const int& lda, double w[],
                 int& info) {
  // check optimal size of lwork
  int lwork = -1;
  std::complex<double>* work  = new std::complex<double>[1];
  double* rwork = new double[std::max(1, 3*n-2)];
  LOOPER_FCALL(zheev,ZHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  lwork = static_cast<int>(std::real(work[0]));
  delete [] work;

  // do real work
  work  = new std::complex<double>[lwork];
  LOOPER_FCALL(zheev,ZHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  delete [] rwork;
  delete [] work;
}

inline void gels(const char& trans, const int& m, const int& n,
                 const int& nrhs,
                 float a[], const int& lda,
                 float b[], const int& ldb,
                 int& info) {
  int lwork = std::min(m,n) + std::max(std::max(1,m), std::max(n,nrhs)) * 64;
  float* work  = new float[lwork];
  LOOPER_FCALL(sgels,SGELS)(trans, m, n, nrhs, a, lda, b, ldb,
                            work, lwork,info);
  delete [] work;
}

inline void gels(const char& trans, const int& m, const int& n,
                 const int& nrhs,
                 double a[], const int& lda,
                 double b[], const int& ldb,
                 int& info) {
  int lwork = std::min(m,n) + std::max(std::max(1,m), std::max(n,nrhs)) * 64;
  double* work  = new double[lwork];
  LOOPER_FCALL(dgels,DGELS)(trans, m, n, nrhs, a, lda, b, ldb,
                            work, lwork,info);
  delete [] work;
}

} // namespace lapack_dispatch

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
  lapack_dispatch::syev(jobz, uplo,
                        vector_helper<Vector>::size(w),
                        matrix_helper<Matrix>::begin_ptr(a),
                        vector_helper<Vector>::size(w),
                        vector_helper<Vector>::begin_ptr(w),
                        info);
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
      lapack_dispatch::gels(trans, m, n, 1,
                            matrix_helper<Matrix>::begin_ptr(a), m,
                            vector_helper<Vector>::begin_ptr(b), m,
                            info);
      std::copy(b.begin(), b.begin() + n, x.begin());
    } else {
      std::copy(b.begin(), b.begin() + m, x.begin());
      lapack_dispatch::gels(trans, m, n, 1,
                            matrix_helper<Matrix>::begin_ptr(a), m,
                            vector_helper<Vector>::begin_ptr(x), n,
                            info);
    }
  } else {
    trans = 'T';
    if (m >= n) {
      lapack_dispatch::gels(trans, n, m, 1,
                            matrix_helper<Matrix>::begin_ptr(a), n,
                            vector_helper<Vector>::begin_ptr(b), m,
                            info);
      std::copy(b.begin(), b.begin() + n, x.begin());
    } else {
      std::copy(b.begin(), b.begin() + m, x.begin());
      lapack_dispatch::gels(trans, n, m, 1,
                            matrix_helper<Matrix>::begin_ptr(a), n,
                            vector_helper<Vector>::begin_ptr(x), n,
                            info);
    }
  }
  if (info != 0) throw std::runtime_error("failed in gels");
}

} // end namespace looper

#endif // LOOPER_LAPACK_H
