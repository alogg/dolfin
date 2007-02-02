// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-30
// Last changed: 2007-01-30
//
// Testing experimental code for ordering of mesh entities.

#include <dolfin.h>
#include "test.h"

using namespace dolfin;

int main()
{
  UnitSquare mesh(2, 2);
//  UnitCube mesh(1, 1, 1);

  mesh.init();

  Function f;

  test::LinearForm L(f);
  FEM::disp(mesh, L.test());

  MeshOrdering::order(mesh);
  FEM::disp(mesh, L.test());
}
