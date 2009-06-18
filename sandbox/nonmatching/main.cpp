// Copyright (C) 2009 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-06-17
// Last changed: 

//
// This program demonstrates the interpolation of functions on non-matching 
// meshes.
//

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
  UnitSquare mesh1(64, 64);

  P1::FunctionSpace V0(mesh0);
  P1::FunctionSpace V1(mesh1);

  MyFunction f0(V0);
  //f0.interpolate();

  Function f1(V1);
  f1.interpolate(f0);

  //Vector x; 
  //V1.interpolate(x, my_function0, "non-matching");
  //Function my_function1(V1, x);
  //my_function1.interpolate(x, V1);

  File file0("test_0_.pvd");
  File file1("test_1_.pvd");
  file0 << f0;
  file1 << f1;

  plot(f0);
  plot(f1);

  return 0;
}
