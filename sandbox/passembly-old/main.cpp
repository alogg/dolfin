// Copyright (C) 2006-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-12-01
// Last changed: 2007-05-29
//
// This file is used for testing parallel assembly

#include "partition.h"
#include "PdofMap.h"
#include "Poisson2D.h"
#include <dolfin.h>

using namespace dolfin;
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  cout << "Starting parallel assemble/solve test." << endl;    

  // Initialise PETSc  
  //PETScManager::init(argc, argv);

  // Initializes MPI PETSc
  dolfin::SubSystemsManager::initMPI();
  dolfin::SubSystemsManager::initPETSc();
  // Get number of processes
  int num_processes_int;
  MPI_Comm_size(PETSC_COMM_WORLD, &num_processes_int);
  unsigned int num_processes = num_processes_int;

  // Get this process number
  int process_int;
  MPI_Comm_rank(PETSC_COMM_WORLD, &process_int);
  unsigned int this_process = process_int;

  // Create mesh and partion
  UnitSquare mesh(50, 50);
  MeshFunction<dolfin::uint> cell_partition_function;
  mesh.partition(num_processes, cell_partition_function);

  /*
  if(this_process == 0)
    plot(cell_partition_function);
    */

  // Create linear and bilinear forms
  Function f(mesh, 1.0);
  Poisson2DBilinearForm a; 
  Poisson2DLinearForm L(f); 

  // Prepare for assembly
  DofMapSet dof_map_set_a;
  DofMapSet dof_map_set_L;
  dof_map_set_a.update(a.form(), mesh);
  dof_map_set_L.update(L.form(), mesh);
  UFC ufc_a(a.form(), mesh, dof_map_set_a);
  UFC ufc_L(L.form(), mesh, dof_map_set_L);

  // Initialize global parallel tensor
  dolfin::uint A_size0 = ufc_a.global_dimensions[0];
  dolfin::uint A_size1 = ufc_a.global_dimensions[1];
  dolfin::uint L_size  = ufc_L.global_dimensions[0];
  
  // Create dof mapping for parallel assembly
  PdofMap pdof_map(A_size0, num_processes);
  pdof_map.create(num_processes, this_process, mesh, ufc_a, cell_partition_function);

  // Get dimension of local (this process) matrix
  unsigned int process_dimension[2];
  process_dimension[0] = pdof_map.process_dimensions(0, this_process);
  process_dimension[1] = pdof_map.process_dimensions(1, this_process);
  
  cout << "Dimensions on process " << this_process << endl;
  cout << process_dimension[0] << " " << process_dimension[1] << endl;

  // Create PETSc parallel vectors
  Vec b, x;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, L_size, &b);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, L_size, &x);
  VecZeroEntries(b);

  // Create PETSc parallel matrix with a guess for number of non-zeroes (10 in thise case)
  Mat A;
  MatCreateMPIAIJ(PETSC_COMM_WORLD, process_dimension[0], process_dimension[1], 
                          A_size0, A_size1, 10, PETSC_NULL, 10, PETSC_NULL, &A); 
  MatZeroEntries(A);

  // Create array for modified dofs
  unsigned int* dofs = new unsigned int[ufc_a.local_dimensions[0]];
  for (unsigned int i = 0; i < ufc_a.local_dimensions[0]; i++)
      dofs[i] = 0;

  tic();
  // Assemble over cells
  cout << "Starting assembly on processor " << this_process << endl;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    if (cell_partition_function.get(*cell) != static_cast<unsigned int>(this_process))
      continue;

    // Update to current cell
    ufc_a.update(*cell);

    // Interpolate coefficients on cell
    //for (dolfin::uint i = 0; i < coefficients.size(); i++)
    //  coefficients[i]->interpolate(ufc.w[i], ufc_a.cell, *ufc_a.coefficient_elements[i], *cell);
    
    // Tabulate unmodified dofs for each dimension
    for (dolfin::uint i = 0; i < ufc_a.form.rank(); i++)
      dof_map_set_a[i].tabulate_dofs(ufc_a.dofs[i], ufc_a.cell);
      //dof_mapa[i]->tabulate_dofs(ufc_a.dofs[i], ufc_a.mesh, ufc_a.cell);
      //ufc_a.dof_maps[i]->tabulate_dofs(ufc_a.dofs[i], ufc_a.mesh, ufc_a.cell);

    // Compute parallel dof map
    pdof_map.update(dofs, ufc_a.dofs[0], ufc_a.local_dimensions[0]);

    // Tabulate cell tensor
    ufc_a.cell_integrals[0]->tabulate_tensor(ufc_a.A, ufc_a.w, ufc_a.cell);

    // Add entries to global tensor
    MatSetValues(A, ufc_a.local_dimensions[0], reinterpret_cast<int*>(dofs), 
                  ufc_a.local_dimensions[0], reinterpret_cast<int*>(dofs), ufc_a.A, ADD_VALUES);
    //A.add(ufc_a.A, ufc.local_dimensions, dofs);
  }
  cout << "Finished assembly " << this_process << "  " << toc() << endl;

  delete [] dofs;

  cout << "Starting finalise assmebly " << this_process << endl;
  tic();
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  double time2 = toc();
  cout << "Finished finalise assembly " << this_process << "  " << time2 << endl;

  // View matrix
//  MatView(A, PETSC_VIEWER_STDOUT_WORLD);

/*
  tic();
  Matrix B;
  Assembler assembler;
  assembler.assemble(B, a, mesh);
  double time1 = toc();
  cout << "Standard assemble time " << time1 << endl;
*/

/*
  int i = 0;
  std::set<int> remapped_dof;
  std::pair<std::set<int>::const_iterator, bool> set_return;
  for(dolfin::uint process = 0; process < num_processes; ++process)
    for(CellIterator cell(mesh); !cell.end(); ++cell)
      if(cell_partition_function.get(*cell) == static_cast<unsigned int>(process) )
        for(VertexIterator vertex(cell); !vertex.end(); ++vertex)
        {  
          set_return = remapped_dof.insert(vertex->index());
          if(set_return.second) // If insertion is successful, renumber. Otherwise dof has already been renumbered
            new_map[ vertex->index() ] = i++;
        }
  remapped_dof.clear();

  // Apply some boundary conditions so that the system can be solved
  // Just apply homogeneous Dirichlet bc to the first three vertices

  IS is = 0;
  int nrows = 3;
//  int rows[3] = {0, 1, 2};
  int rows[3] = {new_map[0], new_map[1], new_map[2]};
  ISCreateGeneral(PETSC_COMM_WORLD, nrows, rows, &is);
  PetscScalar one = 1.0;
  MatZeroRowsIS(A, is, one);
  ISDestroy(is);

  double bc_values[3] = {0, 0, 0};
  VecSetValues(b, nrows, rows, bc_values, INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  // Solve linear system system
  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);

  cout << "Starting solve " << process << endl;
  KSPSolve(ksp, b, x);
  cout << "Finished solve " << process << endl;

  // Print parallel  results
//  cout << "Parallel RHS vector " << endl;
//  VecView(b, PETSC_VIEWER_STDOUT_WORLD);
//  cout << "Parallel matrix " << endl;
//  MatView(A, PETSC_VIEWER_STDOUT_WORLD);


  if(process == 0)
    cout << "*** Parallel solution vector " << endl;
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  if(process == 0)
  {
    set("output destination", "silent");
    PETScMatrix Aref;
    PETScVector bref, xref;
    FEM::assemble(a, L, Aref, bref, mesh); 
    Aref.ident(rows, nrows);
    bref(0) = 0.0; bref(1) = 0.0; bref(2) = 0.0;

    KrylovSolver solver;
    solver.solve(Aref, xref, bref);
    set("output destination", "terminal");

  // Print reference results
//    cout << "Single process reference vector " << endl;
//    VecView(bref.vec(), PETSC_VIEWER_STDOUT_SELF);
//    cout << "Single process reference matrix " << endl;
//    MatView(Aref.mat(), PETSC_VIEWER_STDOUT_SELF);
    cout << "*** Single process solution vector " << endl;
    VecView(xref.vec(), PETSC_VIEWER_STDOUT_SELF);

    cout << "Finsished parallel assemble/solve test." << endl;    
  }
*/
  return 0;
}
