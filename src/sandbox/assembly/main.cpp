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
#include <dolfin/UFC.h>

using namespace dolfin;

class DirichletBC : public BoundaryCondition
{
  void eval(BoundaryValue& value, const Point& p, unsigned int i)
  {
  }
};

int main()
{
  UnitSquare mesh(2, 2);

  // Old assembly
  cout << "---------- Old assembly, DOLFIN Matrix ----------" << endl;
  Matrix A;
  PoissonOld::BilinearForm a;
  dolfin_log(false);
  tic();
  FEM::assemble(a, A, mesh);
  real t0 = toc();
  dolfin_log(true);
  //A.disp();

  // New assembly, DOLFIN matrix
  cout << "---------- New assembly, DOLFIN Matrix ----------" << endl;
  Matrix B;
  Poisson b;
  tic();
  assemble(B, b, mesh);
  real t1 = toc();
  //B.disp();

  // New assembly
  cout << "---------- New assembly, AssemblyMatrix ----------" << endl;
  AssemblyMatrix C;
  Poisson c;
  tic();
  assemble(C, c, mesh);
  real t2 = toc();
  //C.disp();

  cout << "Old assembly, DOLFIN Matrix:  " << t0 << endl;
  cout << "New assembly, DOLFIN Matrix:  " << t1 << endl;
  cout << "New assembly, AssemblyMatrix: " << t2 << endl;

  cout << "---------- Sparsity pattern ----------" << endl;
  DofMaps dof_maps;
  dof_maps.update(b, mesh);
  SparsityPattern sparsity_pattern; 
 
  tic();
  dof_maps.sparsityPattern(sparsity_pattern);
  real t3  = toc();
  cout << "Sparsity pattern: " << t3 << endl;

  uBlasSparseMatrix ublas_matrix;
  cout << "------ Initialise uBlas matrix with sparsity pattern ------" << endl;
  tic();
  ublas_matrix.init(sparsity_pattern);
  real t4  = toc();
  cout << "Matrix initialisation with sparsity pattern: " << t4 << endl;

  DofMaps dof_maps2;
  dof_maps2.update(c, mesh);
  UFC ufc(c, mesh, dof_maps2);

  DirichletBC bc;
  bc.applyBC(B, mesh, **(ufc.finite_elements), **(ufc.dof_maps));
  
  File file("mesh.pvd");
  file << mesh;


}
