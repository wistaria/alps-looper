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

/* $Id: lapack.h 693 2004-03-16 15:48:04Z wistaria $ */

#ifndef LOOPER_LAPACK_H
#define LOOPER_LAPACK_H

#include <looper/fortran.h>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <stdexcept>
#include <vector>

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

// Simple Driver Routines for Eigenvalue Problems for
// Symmetric/Hermitian Matrixes

// SSYEV: compute all eigenvalues and, optionally, eigenvectors of a
// real symmetric matrix (single precision, simple driver)
void LOOPER_FCALL(ssyev,SSYEV)(const char& jobz, const char& uplo, const int& n,
                             float a[], const int& lda,
                             float w[],
                             float work[],
                             const int& lwork, int& info);
  
// DSYEV: compute all eigenvalues and, optionally, eigenvectors of a
// real symmetric matrix (double precision, simple driver)
void LOOPER_FCALL(dsyev,DSYEV)(const char& jobz, const char& uplo, const int& n,
                             double a[], const int& lda,
                             double w[],
                             double work[],
                             const int& lwork, int& info);

// CHEEV: compute all eigenvalues and, optionally, eigenvectors of a
// complex Hermitian matrix (single precision, simple driver)
void LOOPER_FCALL(cheev,CHEEV)(const char& jobz, const char& uplo, const int& n,
                             std::complex<float> a[], const int& lda,
                             float w[],
                             std::complex<float> work[],
                             const int& lwork,
                             float rwork[],
                             int& info);

// ZHEEV: compute all eigenvalues and, optionally, eigenvectors of a
// complex Hermitian matrix (double precision, simple driver)
void LOOPER_FCALL(zheev,ZHEEV)(const char& jobz, const char& uplo, const int& n,
                             std::complex<double> a[], const int& lda,
                             double w[],
                             std::complex<double> work[],
                             const int& lwork,
                             double rwork[],
                             int& info);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

namespace looper {

namespace lapack_dispatch {

inline void syev(const char & jobz, const char& uplo, const int& n,
                 float a[], const int& lda, float w[],
                 int & info) {
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

inline void syev(const char & jobz, const char& uplo, const int& n,
                 double a[], const int& lda, double w[],
                 int & info) {
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

inline void syev(const char & jobz, const char& uplo, const int& n,
                 std::complex<float> a[], const int& lda, float w[],
                 int & info) {
  // check optimal size of lwork
  int lwork = -1;
  std::complex<float>* work  = new std::complex<float>[1];
  float* rwork = new float[std::max(1, 3*n-2)];
  LOOPER_FCALL(cheev,CHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  lwork = static_cast<int>(real(work[0]));
  delete [] work;

  // do real work
  work  = new std::complex<float>[lwork];
  LOOPER_FCALL(cheev,CHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  delete [] rwork;
  delete [] work;
}

inline void syev(const char & jobz, const char& uplo, const int& n,
                 std::complex<double> a[], const int& lda, double w[],
                 int & info) {
  // check optimal size of lwork
  int lwork = -1;
  std::complex<double>* work  = new std::complex<double>[1];
  double* rwork = new double[std::max(1, 3*n-2)];
  LOOPER_FCALL(zheev,ZHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  lwork = static_cast<int>(real(work[0]));
  delete [] work;

  // do real work
  work  = new std::complex<double>[lwork];
  LOOPER_FCALL(zheev,ZHEEV)(jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
  delete [] rwork;
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

template<class T>
struct matrix_helper<boost::numeric::ublas::matrix<T> >
{
  typedef T                                value_type;
  typedef boost::numeric::ublas::matrix<T> matrix_type;
  static value_type * begin_ptr(matrix_type& m) { return &m(0, 0); }
};

}

template <class Matrix, class Vector>
inline void diagonalize(Matrix& a, Vector& w, bool need_eigenvectors = true)
{
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

} // end namespace looper

#endif // LOOPER_LAPACK_H
