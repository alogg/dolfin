// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Simple test program for solving a singular system
// with GMRES + AMG (PETSc + Hypre BoomerAMG)
//
// Singular solver converges much better

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

    real eval(const real* x) const
    {
      real dx = x[0] - 0.5;
      real dy = x[1] - 0.5;
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

  // Compute constant, c = int_\Omega b and subtract it from b  
  real c = 0.0;
  for (unsigned int i=0; i<b.size(); i++)
    c += b[i];
  c /= static_cast<real>(b.size());
  Vector bhat(b.size());
  for (unsigned int i=0; i<b.size(); i++)
    bhat.setitem(i, b[i] - c);

  // Solve modified linear system using ordinary linear solver
  LinearSolver s2(gmres, amg);
  s2.solve(A, x, bhat);
  cout << "Residual: " << residual(A, x, b) << endl;

  return 0;
}
