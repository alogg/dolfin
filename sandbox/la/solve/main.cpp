// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2007-08-20
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 500*exp(-((x-0.5)^2 + (y-0.5)^2)/0.02)
//
// and boundary conditions given by
//
//     u(x, y)     = 0               for x = 0
//     du/dn(x, y) = 25 sin(5 pi y)  for x = 1
//     du/dn(x, y) = 0               otherwise

#include <dolfin.h>
#include "Poisson.h"
  
using namespace dolfin;

int main()
{
  // Source term
  class Source : public Function
  {
  public:
    
    Source(Mesh& mesh) : Function(mesh) {}

    double eval(const double* x) const
    {
      double dx = x[0] - 0.5;
      double dy = x[1] - 0.5;
      return 500.0*exp(-(dx*dx + dy*dy)/0.02);
    }

  };

  // Neumann boundary condition
  class Flux : public Function
  {
  public:

    Flux(Mesh& mesh) : Function(mesh) {}

    double eval(const double* x) const
    {
      if (x[0] > DOLFIN_EPS)
        return 25.0*sin(5.0*DOLFIN_PI*x[1]);
      else
        return 0.0;
    }

  };

  // Sub domain for Dirichlet boundary condition
  class DirichletBoundary : public SubDomain
  {
    bool inside(const double* x, bool on_boundary) const
    {
      return x[0] < DOLFIN_EPS && on_boundary;
    }
  };

  // Create mesh
  UnitSquare mesh(512, 512);

  // Create functions
  Source f(mesh);
  Flux g(mesh);

  // Create boundary condition
  Function u0(mesh, 0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);
  
  // Define PDE
  PoissonBilinearForm a;
  PoissonLinearForm L(f, g);
  LinearPDE pde(a, L, mesh, bc);

  // Assemble matrix and vector
  Matrix A;
  Vector b;
  assemble(A, a, mesh);
  assemble(b, L, mesh);
  bc.apply(A, b, a);
  
  // Solve linear system
  Vector x;
  Function u(mesh, x, a);

  // First option
  tic();
  for (int i = 0; i < 10; i++)
    solve(A, u.vector(), b, gmres, amg);
  message("--- Calling solve repeatedly: %g seconds", toc());

  // Second option
  tic();
  KrylovSolver krylov_solver(gmres, amg);
  for (int i = 0; i < 10; i++)
    krylov_solver.solve(A, u.vector(), b);
  message("--- Reusing KrylovSolver: %g seconds", toc());

  // Third option
  tic();
  LinearSolver linear_solver(gmres, amg);
  for (int i = 0; i < 10; i++)
    linear_solver.solve(A, u.vector(), b);
  message("--- Reusing LinearSolver: %g seconds", toc());
  
  // Plot solution
  plot(u);

  // Save solution to file
  File file("poisson.pvd");
  file << u;

  return 0;
}
