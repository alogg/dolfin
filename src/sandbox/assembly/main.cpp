// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-17
// Last changed: 2007-03-09
//
// Testing experimental code for new assembly.

#include <dolfin.h>

#include "PoissonOld.h"
#include "Poisson.h"

using namespace dolfin;

int main()
{
  UnitSquare mesh(256, 256);

  // Old assembly
  cout << "---------- Old assembly ----------" << endl;
  Matrix A;
  PoissonOld::BilinearForm a;
  dolfin_log(false);
  tic();
  FEM::assemble(a, A, mesh);
  real t0 = toc();
  dolfin_log(true);
  //A.disp();

  // New assembly
  cout << "---------- New assembly ----------" << endl;
  AssemblyMatrix B;
  Poisson b;
  tic();
  assemble(B, b, mesh);
  real t1 = toc();
  //B.disp();

  cout << "Old assembly: " << t0 << endl;
  cout << "New assembly: " << t1 << endl;

  cout << "---------- Sparsity pattern ----------" << endl;
  DofMaps dof_maps;
  dof_maps.update(b, mesh);
  SparsityPattern sparsity_pattern; 
 
  tic();
  dof_maps.sparsityPattern(sparsity_pattern);
  real t2  = toc();

  cout << "Sparsity pattern: " << t2 << endl;

  dolfin::uint* nzrow = new dolfin::uint[sparsity_pattern.size(0)];
  sparsity_pattern.numNonZeroPerRow(nzrow);
  delete [] nzrow;  

}
