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

    virtual double integrate(const int n);
};

class IntegrateSimpson : public Integrate
{
  private:
    double a, b;
    double (*f)(double);

  public:
    IntegrateSimpson(const double a, const double b, double (*f)(double))
     : a(a), b(b), f(f) {}

    virtual double integrate(const int n);
};

class IntegrateRectangle : public Integrate
{
  private:
    double a, b;
    double (*f)(double);

  public:
    IntegrateRectangle(const double a, const double b, double (*f)(double))
     : a(a), b(b), f(f) {}

    virtual double integrate(const int n);
};



double integrate(const double a, const double b, const double tol, std::string method, double (*f)(double));

void print_instructions();

