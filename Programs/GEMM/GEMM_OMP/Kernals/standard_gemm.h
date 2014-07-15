//
// author: Ed Valeev (eduard@valeyev.net)
// date  : July 9, 2014
// the use of this software is permitted under the conditions GNU General Public License (GPL) version 2
//

#ifndef __s2i2_core_gemmkernel_h_DEFINED
#define __s2i2_core_gemmkernel_h_DEFINED

// standard C++ headers
#include <cstddef>

void dgemm(const double* a, const double* b, double* c, size_t n);

void dgemm_blas(const double* a, const double* b, double* c, const size_t n);

#ifdef HAVE_EIGEN
void dgemm_eigen(const double* a, const double* b, double* c, size_t n);
#endif

double matrix_difference(const int N, const double* a, const double* bench);

#endif // __s2i2_core_gemmkernel_h_DEFINED
