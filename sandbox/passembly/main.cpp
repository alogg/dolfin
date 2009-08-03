// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2008-12-26
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 500*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)
//
// and boundary conditions given by
//
//     u(x, y) = 0 for x = 0 or x = 1

#include <dolfin.h>
#include "PoissonP1.h"
#include "PoissonP2.h"

using namespace dolfin;

namespace Poisson = PoissonP1;
//namespace Poisson = PoissonP2;

// Source term
class Source : public Function
{
  void eval(double* values, const double* x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 500.0*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const double* x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

int main()
{
  // Create mesh and function space
  //Mesh mesh("unitsquare.xml.gz");
  Mesh mesh("unitsquare_reallysmall.xml.gz");
  info(mesh.data());

  Poisson::FunctionSpace V(mesh);
  info(mesh.data());

  // Define boundary condition
  Constant u0(0.0);
  DirichletBoundary boundary;
  DirichletBC bc(V, u0, boundary);

  // Define variational problem
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f;
  L.f = f;

  // Compute solution
  VariationalProblem problem(a, L, bc);
  problem.parameters("linear_solver") = "iterative"; // Avoid direct solver for now, seems to break
  Function u;
  problem.solve(u);
  cout << "Process num " << dolfin::MPI::process_number() << endl;
  u.vector().disp();

  double norm = u.vector().norm("l2");
  if (dolfin::MPI::process_number() == 0)
    cout << "Norm of solution vector: " << norm << endl;

  // Plot solution
  //plot(u);

  // Save solution in VTK format
  File file("result/output.pvd");
  file << u;

  return 0;
}
