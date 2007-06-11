// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-01-30
// Last changed: 2007-02-27
//
// Testing experimental code for ordering of mesh entities.

#include <dolfin.h>
#include "Poisson3D_2.h"

using namespace dolfin;

int main()
{
// Right-hand side, 3D
class Source3D : public Function
{
public:
  
  Source3D(Mesh& mesh) : Function(mesh) {}

  real eval(const real* x) const
  {
    return 3.0*DOLFIN_PI*DOLFIN_PI*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])*sin(DOLFIN_PI*x[2]);
  }

};

//  UnitSquare mesh(5, 5);
  UnitCube mesh(1, 1, 1);
//  mesh.order();
//  mesh.init();
//  mesh.init(0);
//  mesh.init(1);
//  mesh.init(2);
//  mesh.init(3);
//  mesh.init(4);

//  mesh.disp();
//  mesh.order();
  Source3D f(mesh);

  Matrix A;
  Vector b;

  Poisson3D_2BilinearForm a;
  Poisson3D_2LinearForm L(f);

//  assemble(A, a, mesh);
  assemble(b, L, mesh);

  mesh.disp();

//  File file("Adev.xml");
//  file << A;

  File file_b("bdev.xml");
  file_b << b;
}




