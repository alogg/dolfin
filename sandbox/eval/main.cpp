// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-03-11
// Last changed: 2007-03-11
//
// Testing evaluation at arbitrary points

#include <dolfin.h>

using namespace dolfin;

class F : public Function
{
public:
  
  F(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = sin(3.0*x[0])*sin(3.0*x[1]);
  }
  
};

int main()
{
  UnitSquare mesh(10, 10);
  real x[2] = {0.3, 0.3};
  real v[1] = {0.0};

  // A user-defined function
  F f(mesh);

  f.eval(v, x);
  message("f(x) = %g", v[0]);
}
