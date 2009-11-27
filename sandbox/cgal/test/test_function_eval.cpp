// Copyright (C) 2009 Andre Massing
// Licensed under the GNU LGPL Version 2.1.
// 
// First added:  2009-11-07
// Last changed: 2009-11-10
//
// Added to generate wiht GTS a list of reference points and function values
// to check function evaluation with CGAL against it.

#include <dolfin.h>
#include "P1.h"
#include "P3.h"

using namespace dolfin;
using dolfin::uint;

#include "P3.h"

#ifdef HAS_GTS

class F : public Expression
{
public:

  F() : Expression(3) {}

  void eval(double* values, const double* x) const
  {
    values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2]);
  }

};

#elif  HAS_GCAL

class F : public Expression
{
public:

  F() : Expression(3) {}

  void eval(double* values, const double* x) const
  {
    values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2]);
  }

};

#endif

int main()
{
  not_working_in_parallel("This demo");
  UnitCube mesh(20, 20, 20);

  // Create function spaces
  P3::FunctionSpace V0(mesh);

  // Create functions
  Function f0(V0);
  
  F f;
  f0.interpolate(f);
  
  double value = 0.0;
  double X[3]  = {0.0, 0.0, 0.0};
  srand(1);
  uint num_points  = 100000;
  for (uint i = 1; i <= num_points; ++i)
  {
    X[0] = std::rand()/static_cast<double>(RAND_MAX);
    X[1] = std::rand()/static_cast<double>(RAND_MAX);
    X[2] = std::rand()/static_cast<double>(RAND_MAX);
    f.eval(&value, X);
    info("x = %.12e\ty = %.12e\tz = %.12e\tf(x) = %.12e",X[0],X[1],X[2],value);
  }

return 0;
}
