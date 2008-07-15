// Copyright (c) 2008 Nuno Lopes
// Licensed under the GNU LGPL Version 2.1
//
// First added:  2008-06-20

#include <dolfin.h>

using namespace dolfin;

int main()
{
  // Illustrates problem with smoothing near non-convex corners of the mesh
  UnitCircle mesh(40, UnitCircle::left, UnitCircle::sumn);
  for (int i = 0; i < 5; i++)
    mesh.smooth();
  plot(mesh);

  return 0;
}
