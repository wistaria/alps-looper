// dpQLM: lapack.h

#ifndef LOOPER_LAPACK_H
#define LOOPER_LAPACK_H

#include <looper/config.h>
#include <algorithm>
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


template <class Matrix, class Vector>
inline void diagonalize(Matrix& a, Vector& w, bool need_eigenvectors = true)
{
  char jobz;
  if (need_eigenvectors) {
    jobz = 'V';
  } else {
    jobz = 'N';
  }
  char uplo = 'U';
  int info;

  // call dispatcher
  lapack_dispatch::syev(jobz, uplo, w.size(), &a(0,0), w.size(), &w[0], info);
  if (info != 0) throw std::runtime_error("failed in syev");
}

} // end namespace looper

#endif // LOOPER_LAPACK_H
