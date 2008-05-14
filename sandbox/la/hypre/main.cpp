// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Simple test program for solving a singular system
// with GMRES + AMG (PETSc + Hypre BoomerAMG)

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
  UnitSquare mesh(32, 32);

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

  // Solve linear sysem
  solve(A, x, b, gmres, amg);

  return 0;
}
