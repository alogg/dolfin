// Place for random tests

#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

class MyFunction : public Function
{
public:

  MyFunction(const FunctionSpace& V) : Function(V) {};

  void eval(double* values, const double* x) const
  {
    message("Calling eval");
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 500.0*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

class FunctionContainer
{
public:
  FunctionContainer(const FunctionSpace& V)
  {
    Function g(V);
    g.vector();
    _f = g;
  };
	
  const Function& get_function()
  {
    return _f;
  };
protected:
  Function _f;
};


int main()
{  
  UnitSquare mesh(2, 2);
  PoissonFunctionSpace V(mesh);
  MyFunction f(V);
  Vector x;
  
  message("Interpolating to another vector");
  f.interpolate(x, f.function_space());
  x.disp();

  message("Interpolating to the function vector");
  f.interpolate(f.vector(), f.function_space());
  f.vector().disp();
  
  message("Interpolating using initialising by an external function");
  MyFunction f_(f);
  f.interpolate(f_.vector(), f.function_space());
  f.vector().disp();
  
  message("Returning Function by reference");
  FunctionContainer fc(V);
  const Function& f2 = fc.get_function();

  message("dim = %d", f2.function_space().dim());
}

