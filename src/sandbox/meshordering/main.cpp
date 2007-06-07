// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-01-30
// Last changed: 2007-02-27
//
// Testing experimental code for ordering of mesh entities.

#include <dolfin.h>

using namespace dolfin;

int main()
{
//  UnitSquare mesh(2, 2);
  UnitCube mesh(2, 2, 2);
  mesh.order();
//  mesh.init();
//  mesh.init(0);
  mesh.init(1);
//  mesh.init(2);
  mesh.init(3);
//  mesh.init(4);

//  mesh.disp();
//  mesh.order();
  mesh.disp();
}
