//
// author: Daniel Smith (dsmith@auburn.edu)
// date  : July 14, 2014
// the use of this software is permitted under the conditions GNU General Public License (GPL) version 2
//

//#include <cmath>
#include <omp.h>
#include <iostream>

#include "blocked_gemm.h"


// The inner gemm component
void ikj_simple_gemm(double* a, double* b, double* c, size_t n,
               int starti, int endi,
               int startj, int endj,
               int startk, int endk){

  for (int i=starti; i < endi; i++){
    for (int k=startk; k < endk; k++){
      const int in = i * n;
      const int kn = k * n;
      const double aik = a[in + k];

#pragma simd
      for (int j=startj; j < endj; j++){
        c[in + j] += aik * b[kn + j];
      }
    }
  }
}

// The inner gemm component
void kij_simple_gemm(double* a, double* b, double* c, size_t n,
               int starti, int endi,
               int startj, int endj,
               int startk, int endk){

  for (int k=startk; k < endk; k++){
    for (int i=starti; i < endi; i++){
      const int kn = k*n;
      const int in = i*n;
      const double aik = a[in + k];

#pragma simd
      for (int j=startj; j < endj; j++){
        c[in + j] += aik * b[kn + j];
      }
    }
  }
}

// The inner gemm component
void ikj_unroll_gemm(double* a, double* b, double* c, size_t n,
               int starti, int endi,
               int startj, int endj,
               int startk, int endk){

  for (int i=starti; i < endi; i++){
    for (int k=startk; k < endk; k+=4){
      const int in = i*n;
      const int kn0 = k*n;
      const int kn1 = kn0 + n;
      const int kn2 = kn1 + n;
      const int kn3 = kn2 + n;

      const double aik0 = a[in + k];
      const double aik1 = a[in + k + 1];
      const double aik2 = a[in + k + 2];
      const double aik3 = a[in + k + 3];

#pragma simd
#pragma prefetch c:0:4
      for (int j=startj; j < endj; j+=4){
        c[in + j] += aik0 * b[kn0 + j] + aik1 * b[kn1 + j]
                   + aik2 * b[kn2 + j] + aik3 * b[kn3 + j];
        c[in + j+1] += aik0 * b[kn0 + j+1] + aik1 * b[kn1 + j+1]
                    + aik2 * b[kn2 + j+1] + aik3 * b[kn3 + j+1];
        c[in + j+2] += aik0 * b[kn0 + j+2] + aik1 * b[kn1 + j+2]
                    + aik2 * b[kn2 + j+2] + aik3 * b[kn3 + j+2];
        c[in + j+3] += aik0 * b[kn0 + j+3] + aik1 * b[kn1 + j+3]
                    + aik2 * b[kn2 + j+3] + aik3 * b[kn3 + j+3];
      }
    }
  }

}


// The inner gemm component
void kij_unroll_gemm(double* a, double* b, double* c, size_t n,
               int starti, int endi,
               int startj, int endj,
               int startk, int endk){

  for (int k=startk; k < endk; k++){
    for (int i=starti; i < endi; i+=4){
      const int kn = k*n;

      const int in0 = i*n;
      const int in1 = in0+n;
      const int in2 = in1+n;
      const int in3 = in2+n;

      const double aik0 = a[in0 + k];
      const double aik1 = a[in1 + k];
      const double aik2 = a[in2 + k];
      const double aik3 = a[in3 + k];

#pragma simd
#pragma prefetch c:1:4
      for (int j=startj; j < endj; j++){
        const double bkj0 = b[kn + j];
        c[in0 + j] += aik0 * bkj0;
        c[in1 + j] += aik1 * bkj0;
        c[in2 + j] += aik2 * bkj0;
        c[in3 + j] += aik3 * bkj0;
      }
    }
  }

}

void dgemm_blocked(double* a, double* b, double* c, size_t n, size_t block_size){

   const size_t i_block_size = block_size*4;
   const size_t j_block_size = block_size*16;
   const size_t k_block_size = block_size;
  // Start for loops

#pragma omp parallel for default(none) shared(a, b, c, n, block_size)
  for(int i = 0; i < n; i+=i_block_size){
    const size_t starti = i;
    size_t endi = starti + i_block_size;
    if (endi>n) endi = n; 

// #pragma omp parallel for default(none) shared(a, b, c, n, block_size, endi, starti)
    for(int j = 0; j < n; j+=j_block_size){
      const size_t startj = j;
      size_t endj = startj + j_block_size;
      if (endj>n) endj = n; 

      for(int k = 0; k < n; k+=k_block_size){
        const size_t startk = k;
        size_t endk = startk + k_block_size;
        if (endk>n) endk = n; 

        // ikj_simple_gemm(a, b, c, n, starti, endi, startj, endj, startk, endk);
        ikj_unroll_gemm(a, b, c, n, starti, endi, startj, endj, startk, endk);

        // kij_simple_gemm(a, b, c, n, starti, endi, startj, endj, startk, endk);
        // kij_unroll_gemm(a, b, c, n, starti, endi, startj, endj, startk, endk);

      }
    }
  }
}




