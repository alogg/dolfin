// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-17
// Last changed: 2007-01-17
//
// Testing experimental code for new assembly.

#include <dolfin.h>

#include "poisson_ufc.h"
#include "AssemblyMatrix.h"

using namespace dolfin;

int main()
{
  AssemblyMatrix A;
  poisson form;
  UnitSquare mesh(5, 5);

  assemble(A, form, mesh);
}
