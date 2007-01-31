// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-30
// Last changed: 2007-01-30
//
// Testing experimental code for ordering of mesh entities.

#include <dolfin.h>

using namespace dolfin;

int main()
{
  UnitSquare mesh(5, 5);

  MeshOrdering::order(mesh);
}
