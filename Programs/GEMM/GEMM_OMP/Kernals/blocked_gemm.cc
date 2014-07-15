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
void simple_inner_gemm(double* a, double* b, double* c, size_t n,
               int starti, int endi,
               int startj, int endj,
               int startk, int endk){

#pragma omp parallel for default(none) shared(a, b, c, n, starti, endi, startj, endj, startk, endk)
    for (int i=starti; i < endi; i++){
#pragma unroll_and_jam(8)
      for (int j=startj; j < endj; j++){
        const int jn = j*n;
        double tmp_c = c[i + jn];
#pragma simd vectorlength(4)
#pragma unroll_and_jam(8)
        for (int k=startk; k < endk; k++){
          tmp_c += a[i + k*n] * b[k + jn];
        }
        c[i + jn] = tmp_c;
      }
    }

}

// The unrolled inner gemm component
// Apparently the intel compiler can do this better than I can with unroll_and_jam
void unrolled_inner_gemm(double* a, double* b, double* c, size_t n,
               int starti, int endi,
               int startj, int endj,
               int startk, int endk){

#pragma omp parallel for default(none) shared(a, b, c, n, starti, endi, startj, endj, startk, endk)
    for (int i=starti; i < endi; i++){
      for (int j=startj; j < endj; j+=4){
        const int jn0 = j*n;
        const int jn1 = jn0 + n;
        const int jn2 = jn1 + n;
        const int jn3 = jn2 + n;

        double tmp_c0 = c[i + jn0];
        double tmp_c1 = c[i + jn1];
        double tmp_c2 = c[i + jn2];
        double tmp_c3 = c[i + jn3];

#pragma simd vectorlength(4)
        for (int k=startk; k < endk; k+=2){
          const int aik1 = a[i + k*n];
          const int aik2 = a[i + k*n + n];

          tmp_c0 += aik1 * b[k + jn0] + aik2 * b[k + jn0 + 1];
          tmp_c1 += aik1 * b[k + jn1] + aik2 * b[k + jn1 + 1];
          tmp_c2 += aik1 * b[k + jn2] + aik2 * b[k + jn2 + 1];
          tmp_c3 += aik1 * b[k + jn3] + aik2 * b[k + jn3 + 1];
        }
        c[i + jn0] = tmp_c0;
        c[i + jn1] = tmp_c1;
        c[i + jn2] = tmp_c2;
        c[i + jn3] = tmp_c3;
      }
    }

}


void dgemm_blocked(double* a, double* b, double* c, size_t n, size_t block_size){

  // Compute number of blocks in each direction
  int tmp_num_blocks;
  if ((n % block_size) == 0){
    tmp_num_blocks = n / block_size;
  }
  else{
    tmp_num_blocks = n / block_size + 1;
  }
  const int num_blocks = tmp_num_blocks;

  // Compute starting and ending arrays
  double* start = new double[num_blocks];
  double* end = new double[num_blocks];

  int block_pos = 0;

// Appears to actually slow down the computation, simply not big enough to matter
//#pragma omp parallel for default(none) shared(start, end, block_size)
  for(int i=0; i < num_blocks; i++){
    start[i] = block_size * i;
    end[i] = block_size * (i + 1);
  }
  
  end[num_blocks-1] = n;

  // Start for loops
  for(int i = 0; i < num_blocks; i++){
    const int starti = start[i];
    const int endi = end[i];
    for(int j = 0; j < num_blocks; j++){
      const int startj = start[j];
      const int endj = end[j];
      for(int k = 0; k < num_blocks; k++){
        const int startk = start[k];
        const int endk = end[k];
        simple_inner_gemm(a, b, c, n, starti, endi, startj, endj, startk, endk);
        // unrolled_inner_gemm(a, b, c, n, starti, endi, startj, endj, startk, endk);
      }
    }
  }
  delete[] start;
  delete[] end; 
}



