//
// author: Ed Valeev (eduard@valeyev.net)
// date  : July 9, 2014
// the use of this software is permitted under the conditions GNU General Public License (GPL) version 2
//
// edited: Daniel Smith (dsmith@auburn.edu)
// data  : July 14, 2014
//


// standard C++ headers
#include <ctime>
#include <numeric>
#include <iostream>
#include <sstream>
#include <cassert>
#include <chrono>

#include "Kernals/standard_gemm.h"
#include "Kernals/blocked_gemm.h"

double profile_dgemm(size_t n,
                     size_t nrepeats,
                     std::string kernel,
                     bool check);

void check_method(double* a, double* b, double* c, size_t n);


int main(int argc, char* argv[]) {

  // validate command-line arguments
  if (argc != 4 && argc != 5) {
    std::cout << "gemm -- benchmarks square matrix multiplication C[i][j] = A[i][k] B[k][j]" << std::endl;
    std::cout << "usage: gemm n nrepeats kernel" << std::endl;
    std::cout << "       n        -- number of elements in vectors x and y" << std::endl;
    std::cout << "       nrepeats -- number of times to run the kernel" << std::endl;
    std::cout << "       kernel   -- kernel type, allowed values:" << std::endl;
    std::cout << "                   plain   = plain ole loops" << std::endl;
    std::cout << "                   blockYY = blocked loops, block size = YY" << std::endl;
    std::cout << "                   blas    = call to BLAS library's dgemm function" << std::endl;
#ifdef HAVE_EIGEN
    std::cout << "                   eigen   = call to Eigen" << std::endl;
#endif
    std::cout << "       check    -- optional: if 1 check the method against BLAS DGEMM" << std::endl;
    std::cout << "                   some algorithms should only be checked when nrepeats=0" << std::endl;
    return 1;
  }

  // If 5 arguements run check
  bool check = 0;
  if (argc == 5){
    check = 1;
  } 

  std::stringstream ss; ss << argv[1] << " " << argv[2] << " " << argv[3];

  size_t n; ss >> n;
  size_t nrepeats; ss >> nrepeats;
  std::string kernel; ss >> kernel;
  assert(kernel == "plain" ||
         kernel == "blas"  ||
         kernel == "eigen" ||
         kernel.find("block") != std::string::npos);

  profile_dgemm(n, nrepeats, kernel, check);

  return 0;
}

void check_method(double* a, double* b, double* c, size_t n){

  // Obtain relative difference between matrices
  size_t nsq = n*n;
  double* bench = new double[nsq];
  std::fill(bench, bench+nsq, 0.0);

  dgemm_blas(a, b, bench, n);
  double reldiff = matrix_difference(n, c, bench);
  std::cout << "The relative difference between the current kernel and C_DGEMM is " << reldiff << "%\n";

  delete[] bench;
} 

double profile_dgemm(size_t n,
                     size_t nrepeats,
                     std::string kernel,
                     bool check) {


  size_t nsq = n*n;
  double* a = new double[nsq];
  double* b = new double[nsq];
  double* c = new double[nsq];
  std::fill(a, a+nsq, 4.0);
  std::fill(b, b+nsq, 2.0);
  std::fill(c, c+nsq, 0.0);

  const auto tstart = std::chrono::system_clock::now();

  if (kernel == "plain"){
    for (int i=0; i<nrepeats; i++){
      dgemm(a, b, c, n);
    }
  }
  else if (kernel == "blas"){
    for (int i=0; i<nrepeats; i++){
      dgemm_blas(a, b, c, n);
    }
  }
#ifdef HAVE_EIGEN
  else if (kernel == "eigen"){
    for (int i=0; i<nrepeats; i++){
      dgemm_eigen(a, b, c, n);
    }
  }
#endif
  else if (kernel.find("block") != std::string::npos) {
    std::stringstream ss; ss << std::string(kernel.begin()+5, kernel.end());
    size_t blocksize; ss >> blocksize;
    for (int i=0; i<nrepeats; i++){
      dgemm_blocked(a, b, c, n, blocksize);
    }
  }
  else {
    std::cerr << "invalid kernel" << std::endl;
    exit(1);
  }

  const auto tstop = std::chrono::system_clock::now();
  const std::chrono::duration<double> time_elapsed = tstop - tstart;

  std::cout << "n = " << n << " nrepeats = " << nrepeats << " kernel = " << kernel
            << " elapsed time = " << time_elapsed.count() << " s"
            << " throughput = " << (2*n + 1) * nsq * nrepeats / (time_elapsed.count() * 1.e9) << " GFLOP/s"<< std::endl;

  // evaluate the "trace" of c ... not reading c may allow
  // the compiler to skip the computation altogether
  const double c_trace = std::accumulate(c, c+nsq, 0.0);

  if (check){
    check_method(a, b, c, n);
  }

  // Cleanup
  delete[] a;
  delete[] b;
  delete[] c;

  return c_trace;
}

