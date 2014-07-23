#include <cstdio>
#include <cmath>
#include <iostream>
#include <chrono>

#include "integrate_kernals.h"

// Our function
double func(double x){
  // double result = x*x;
  double result = x*x*x*x - 5*x*x*x + 3;
  return result;
}


int main(int argc, char** argv) {
  if (argc==1){
    print_instructions();
    return 1;
  }
  std::string method = argv[1];
  
  const double a = 0;
  const double b = 10;
  const double tol = 1E-5;
//  std::string method = "simpson";

  double result = integrate(a, b, tol, method, func);
  return 0;
}

