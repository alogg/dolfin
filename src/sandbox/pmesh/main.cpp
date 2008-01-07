// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-10-30
// Last changed: 2007-12-17
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>
#include "Poisson2D.h"

using namespace dolfin;
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  UnitCube mesh(4, 4, 4);
  MeshFunction<dolfin::uint> partitions;
  mesh.partition(dolfin::MPI::numProcesses(), partitions);

  partitions.disp();

  //Matrix B;
  //Poisson2DBilinearForm a;
  //pAssembler assembler(mesh);
  //assembler.assemble(B, a, true);

  return 0;
}
