// Copyright (C) 2009 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin.h>
#include "P1.h"

using namespace dolfin;

class MyFunction : public Function
{
  public:

  MyFunction(FunctionSpace& V) : Function(V) {}

  void eval(double* values, const double* x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 500.0*exp(-(dx*dx + dy*dy) / 0.02);
  }
};


int main()
{
  UnitSquare mesh0(16, 16);
  UnitCube mesh1(32, 32, 32);

  P1::FunctionSpace V0(mesh0);
  P1::FunctionSpace V1(mesh1);

  MyFunction my_function0(V0);
  my_function0.interpolate();

  Vector x; 
  V1.interpolate(x, my_function0);

  Function my_function1(V1, x);

  //my_function1.interpolate(x, V1);

  plot(my_function0);
  plot(my_function1);

  return 0;
}
