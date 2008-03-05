// Copyright (C) 2007 Anders Logg and Marie Rognes
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-04-20
// Last changed: 2007-04-20
//
// This demo program solves the mixed formulation of
// Poisson's equation:
//
//     q + grad(u) = 0
//          div(q) = f
//
// The corresponding weak (variational problem)
//
//     <v, q> - <div(v), u> = 0       for all v
//               <w, div q> = <w, f>  for all w
//
// is solved using BDM (Brezzi-Douglas-Marini) elements
// of degree 1 (v, q) and DG (discontinuous Galerkin)
// elements of degree 0 for (w, u).

#include <dolfin.h>
#include "MixedPoisson.h"

using namespace dolfin;

int main()
{
  UnitSquare mesh(1, 1);
  Matrix A;
  MixedPoissonBilinearForm a;

  // Assemble matrix
  assemble(A, a, mesh);
  
  // Display matrix values
  A.disp();

  return 0;
}
