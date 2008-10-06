// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Simple test program for solving a singular system
// with GMRES + AMG (PETSc + Hypre BoomerAMG)
//
// This system has a zero eigenvalue with corresponding
// eigenvector v = (1, 1, 1, ...). This gives a compatibility
// condition (corresponding to the integral of the right-hand
// side being zero for the Neumann problem).
//
// To solve the linear system, we must therefore either make the
// system nonsingular by adding a constraint, like zero average for x,
// or modify the right-hand side to satisfy the compatibility
// condition.

#include <dolfin.h>
#include "Poisson.h"
  
using namespace dolfin;

int main(int argc, char* argv[])
{
  // Initialize PETSc with command-line arguments
  SubSystemsManager::initPETSc(argc, argv);

  // Monitor convergence
  dolfin_set("Krylov monitor convergence", true);

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

  // Create mesh
  UnitSquare mesh(2, 2);

  // Create functions
  Source f(mesh);

  // Variational forms
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  // Assemble linear system
  Matrix A;
  Vector x, b;
  assemble(A, a, mesh);
  assemble(b, L, mesh);

  // Solve linear system using ordinary linear solver
  LinearSolver s0(gmres, amg);
  s0.solve(A, x, b);
  cout << "Residual: " << residual(A, x, b) << endl;

  // Solve linear system using special singular solver
  SingularSolver s1(gmres, amg);
  s1.solve(A, x, b);
  cout << "Residual: " << residual(A, x, b) << endl;

  // Solve modified linear system using ordinary linear solver
  LinearSolver s2(gmres, amg);
  normalize(b, average);
  s2.solve(A, x, b);
  cout << "Residual: " << residual(A, x, b) << endl;

  return 0;
}
