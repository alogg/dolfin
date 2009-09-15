// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-15
// Last changed: 2009-09-15
//
// Simple Poisson assembler

#include <cstdlib>
#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main(int argc, char* argv[])
{
  // Parse command-line arguments
  if (argc != 2)
  {
    info("Usage: solve-poisson n");
    return 1;
  }
  int n = atoi(argv[1]);

  // Create mesh and function space
  UnitCube mesh(n, n, n);
  Poisson::FunctionSpace V(mesh);

  // Assemble matrix
  Poisson::BilinearForm a(V, V);
  Matrix A;
  double t = time();
  assemble(A, a);
  t = time() - t;

  // Report timing
  info("TIME: %.5g", t);

  return 0;
}
