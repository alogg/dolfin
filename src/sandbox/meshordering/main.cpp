// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-30
// Last changed: 2007-02-27
//
// Testing experimental code for ordering of mesh entities.

#include <dolfin.h>

using namespace dolfin;

int main()
{
  UnitSquare mesh(1, 1);
  mesh.init();

  mesh.disp();
  
  MeshOrdering::order(mesh);

  mesh.disp();
}
