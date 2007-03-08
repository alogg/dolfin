// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-17
// Last changed: 2007-03-06
//
// Testing experimental code for new assembly.

#include <dolfin.h>

#include "PoissonOld.h"
#include "Poisson.h"
#include "AssemblyMatrix.h"

using namespace dolfin;

int main()
{
  UnitSquare mesh(2, 2);

  // Old assembly
  /*
  cout << "---------- Old assembly ----------" << endl;
  Matrix A;
  PoissonOld::BilinearForm a;
  FEM::assemble(a, A, mesh);
  A.disp();
  */

  // New assembly
  cout << "---------- New assembly ----------" << endl;
  AssemblyMatrix B;
  Poisson b;
  assemble(B, b, mesh);
  B.disp();
}
