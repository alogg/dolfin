// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-17
// Last changed: 2007-03-01
//
// Testing experimental code for new assembly.

#include <dolfin.h>

#include "Poisson.h"
#include "poisson_ufc.h"
#include "AssemblyMatrix.h"

using namespace dolfin;

int main()
{
  UnitSquare mesh(2, 2);

  // Old assembly
  //cout << "---------- Old assembly ----------" << endl;
  //Matrix A;
  //Poisson::BilinearForm a;
  //FEM::assemble(a, A, mesh);
  //A.disp();

  // New assembly
  cout << "---------- New assembly ----------" << endl;
  AssemblyMatrix B;
  poisson b;
  assemble(B, b, mesh);
  B.disp();
}
