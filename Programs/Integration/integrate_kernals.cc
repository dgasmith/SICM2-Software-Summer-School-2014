#include <cstdio>
#include <cmath>
#include <iostream>
#include <chrono>


class Integrate
{
  public:
    virtual double integrate(const int n) = 0;
};

class IntegrateTrapz : public Integrate
{
  private:
    double a, b;
    double (*f)(double);

  public:
    IntegrateTrapz(const double a, const double b, double (*f)(double))
     : a(a), b(b), f(f) {}

    virtual double integrate(const int n)
    {
      const double h = (b-a)/n;
      double sum = 0.0;

      for (int i=1; i<(n-1); i++) {
        sum += f(a + i*h);
      }

      sum += 0.5*(f(b) + f(a));
      sum *= h;
      return sum;
    }
};

class IntegrateSimpson : public Integrate
{
  private:
    double a, b;
    double (*f)(double);

  public:
    IntegrateSimpson(const double a, const double b, double (*f)(double))
     : a(a), b(b), f(f) {}

    virtual double integrate(const int n)
    {
       const double h = (b-a)/n;
       // const double m = n/2;
       double sum = f(a) + f(b);

       for (int i=1; i < n; i+=2){
         sum += 4.0 * f(a + i*h);
       }   
       for (int i=2; i < (n-1); i+=2){
         sum += 2.0 * f(a + i*h);
       }   

       sum *= (h/3.0);
       return sum;
    }
};

class IntegrateRectangle : public Integrate
{
  private:
    double a, b;
    double (*f)(double);

  public:
    IntegrateRectangle(const double a, const double b, double (*f)(double))
     : a(a), b(b), f(f) {}

    virtual double integrate(const int n)
    {
       const double h = (b-a)/n;
       // const double m = n/2;
       double sum = f(a) + f(b);

       for (int i=1; i < (n-1); i++){
         sum += f(a + i*h);
       }  

       sum *= h; 
       return sum;
    }
};


double integrate(const double a, const double b, const double tol, std::string method, double (*f)(double)) {

  Integrate* integrate_scheme = 0;
  if (method=="trapz"){
    integrate_scheme = new IntegrateTrapz(a, b, f);
  }
  else if (method=="simpson"){
    integrate_scheme = new IntegrateSimpson(a, b, f);
  }
  else if (method=="rectangle"){
    integrate_scheme = new IntegrateRectangle(a, b, f);
  }
  else{
    std::cout << "Could not find integration scheme " << method << std::endl;
    //print_instructions();
  }

  int current_N = 10;
  const int max_cycle = 20;

  double* value_array = new double[max_cycle];
  double* error_array = new double[max_cycle];


  for (int i=0; i<max_cycle; i++){
    value_array[i] = integrate_scheme->integrate(current_N);
    if (i==0){
      error_array[i] = 1;
    }
    else{
      error_array[i] = std::abs(value_array[i] - value_array[i-1]);
    }

    std::cout << "N: " << current_N << "  Value: " << value_array[i] << " Error: " << error_array[i] << "\n";

    if ( error_array[i] > tol){
        current_N *= 2;
    }
    else{
      const double ret_val = error_array[i];
      delete[] value_array;
      delete[] error_array;

    return ret_val;
    }     

  } 


}

void print_instructions(){
  std::cout << "Recurisvely integrate a function based off different integration schemes" << std::endl;
  std::cout << "usage: integrate(double a, double b, double tol, string method, double f)" << std::endl;
  std::cout << "       a        -- start of integration range" << std::endl;
  std::cout << "       b        -- end of integration range" << std::endl;
  std::cout << "       tol      -- convergence tolerance between integration recursions" << std::endl;
  std::cout << "       method   -- integration scheme, allowed values:" << std::endl;
  std::cout << "                   trapz   :  trapezoidal rule integration" << std::endl;
  std::cout << "                   simpson :  simpson rule integration" << std::endl;
//  std::cout << "                   boole   :  simpson rule integration" << std::endl;
  std::cout << "       func     -- funciton to integrate" << std::endl;
}

