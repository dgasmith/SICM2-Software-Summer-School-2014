//
// author: Ed Valeev (eduard@valeyev.net)
// date  : July 9, 2014
// the use of this software is permitted under the conditions GNU General Public License (GPL) version 2
// 
// edited: Daniel Smith (dsmith@auburn.edu)
// data : July 14, 2014
// Add OMP parallelization
//


#include <cassert>
//#include <cmath>
#include <omp.h>
#ifdef HAVE_MKL
#  include <mkl_cblas.h>
#else
#  include <cblas.h>
#endif

// Eigen library Core capabilities
#ifdef HAVE_EIGEN
#  include <Eigen/Core>
#endif

#include "standard_gemm.h"

void dgemm(const double* a, const double* b, double* c, size_t n) {

    size_t ij = 0;
    for (int i = 0; i < n; ++i) {
#pragma omp parallel for default(none) shared(a, b, c, n, i) private(ij)
      for (int j = 0; j < n; ++j, ++ij) {

        double v = 0.0;
        size_t ik = i * n;
        size_t kj = j;

#pragma simd vectorlength(4)
#pragma unroll_and_jam(8)
        for (int k = 0; k < n; ++k, ++ik, kj += n) {
          v += a[ik] * b[kj];
        }
        c[ij] = v;
      }
    }
}

void dgemm_blas(const double* a, const double* b, double* c, const size_t n) {

  size_t tmp_n = n; 
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, tmp_n, tmp_n, tmp_n, 1.0, a, tmp_n,
                b, tmp_n, 0.0, c, tmp_n);

}

#ifdef HAVE_EIGEN
void dgemm_eigen(const double* a, const double* b, double* c, size_t n) {

  using namespace Eigen;
  typedef Eigen::Matrix<double,
                        Eigen::Dynamic,
                        Eigen::Dynamic,
                        Eigen::RowMajor> Matrix; // row-major dynamically-sized matrix of double
  Eigen::Map<const Matrix> aa(a, n, n);
  Eigen::Map<const Matrix> bb(b, n, n);
  Eigen::Map<Matrix> cc(c, n, n);
  cc = aa * bb; 
}
#endif

double matrix_difference(const int N, const double* a, const double* bench) {
  double sum = 0.0;
  double divisor;
  const double thresh = 1E-15;
 #pragma omp parallel for default(none) shared(N, a, bench) private(divisor) reduction(+:sum)
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      divisor = bench[i*N+j];
      if (std::abs(divisor) < thresh){
        divisor = thresh;
      }
      // sum += std::abs(a[i*N+j]-bench[i*N+j]);
      sum += std::abs((a[i*N+j]-bench[i*N+j])/divisor);
    }
  }
  sum /= (N*N/100.0);
  return sum;
 }


