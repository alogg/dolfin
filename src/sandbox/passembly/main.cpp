// Copyright (C) 2006 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-12-01
// Last changed: 
//
// This file is used for testing parallel assembly

#include <set>
#include "partition.h"
#include "Poisson2D.h"

using namespace dolfin;
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  cout << "Starting parallel assemble/solve test." << endl;    

  // Initialise PETSc  
  PETScManager::init(argc, argv);

  // Get number of processes
  int num_processes_int;
  MPI_Comm_size(PETSC_COMM_WORLD, &num_processes_int);
  unsigned int num_processes = num_processes_int;

  // Get this process number
  int process_int;
  MPI_Comm_rank(PETSC_COMM_WORLD, &process_int);
  unsigned int this_process = process_int;

  // Create mesh
  UnitSquare mesh(1500, 1500);

  MeshFunction<dolfin::uint> cell_partition_function;
  MeshFunction<dolfin::uint> vertex_partition_function;
  cout << "Creating mesh partition (using METIS)" << endl; 
  tic();
  testMeshPartition(mesh, cell_partition_function, vertex_partition_function,
                    num_processes);
  real time = toc();
  cout << "Finished partitioning mesh " << time << endl;

  // Create linear and bilinear forms
  Function f(mesh, 1.0);
  Poisson2DBilinearForm a; 
  Poisson2DLinearForm L(f); 

  // Prepare for assembly
  DofMaps dof_maps_a;
  DofMaps dof_maps_L;
  dof_maps_a.update(a.form(), mesh);
  dof_maps_L.update(L.form(), mesh);
  UFC ufc_a(a.form(), mesh, dof_maps_a);
  UFC ufc_L(L.form(), mesh, dof_maps_L);

  // Initialize global parallel tensor
  dolfin::uint A_size0 = ufc_a.global_dimensions[0];
  dolfin::uint A_size1 = ufc_a.global_dimensions[1];
  dolfin::uint L_size  = ufc_L.global_dimensions[0];
  
  std::vector<unsigned int> map(A_size0);
  map.clear();
  std::set<unsigned int> set;
  std::pair<std::set<unsigned int>::const_iterator, bool> set_return;

  unsigned int dof_counter = 0;
  std::vector<unsigned int> partition_dof_counter(num_processes);
  partition_dof_counter.clear();

  cout << "Total dofs " << A_size0 << endl;

  for (unsigned int proc = 0; proc < num_processes; ++proc)
  {
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      if (cell_partition_function.get(*cell) != static_cast<unsigned int>(proc))
        continue;

      // Update to current cell
      ufc_a.update(*cell);
  
      // Tabulate dofs for each dimension
      ufc_a.dof_maps[0]->tabulate_dofs(ufc_a.dofs[0], ufc_a.mesh, ufc_a.cell);
      ufc_a.dof_maps[1]->tabulate_dofs(ufc_a.dofs[1], ufc_a.mesh, ufc_a.cell);

      for(unsigned int i=0; i < ufc_a.dof_maps[0]->local_dimension(); ++i)
      {
        set_return = set.insert( (ufc_a.dofs[0])[i] );
        if( set_return.second )
        {
          map[ (ufc_a.dofs[0])[i] ] = dof_counter++;
          partition_dof_counter[proc]++;
        }
      }
    }
  }
  cout << "Number of dofs on process " << this_process << " = " << partition_dof_counter[this_process] << endl;

  // Create PETSc parallel vectors
  Vec b, x;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, L_size, &b);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, L_size, &x);

  // Create PETSc parallel matrix with a guess for number of non-zeroes (10 in thise case)
  Mat A;
  MatCreateMPIAIJ(PETSC_COMM_WORLD, partition_dof_counter[this_process], partition_dof_counter[this_process], 
                  A_size0, A_size1, 10, PETSC_NULL, 10, PETSC_NULL, &A); 

  // Zero matrix and vector
  MatZeroEntries(A);
  VecZeroEntries(b);

  // Assemble over cells
  cout << "Starting assembly on processor " << this_process << endl;
  unsigned int counter = 0;

  // Initialize new dof mapping
/*  int** dofs;
  dofs = new int*[ufc_a.form.rank()];
  for (unsigned int i = 0; i < ufc_a.form.rank(); i++)
  {
    dofs[i] = new int[ufc_a.local_dimensions[i]];
    for (unsigned int j = 0; j < ufc_a.local_dimensions[i]; j++)
      dofs[i][j] = 0;
  }
*/
  int* dofs;
  dofs = new int[ufc_a.local_dimensions[0]];
  for (unsigned int i = 0; i < ufc_a.local_dimensions[0]; i++)
      dofs[i] = 0;

  tic();
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    if (cell_partition_function.get(*cell) != static_cast<unsigned int>(this_process))
      continue;

    // Update to current cell
    ufc_a.update(*cell);

    // Interpolate coefficients on cell
//    for (dolfin::uint i = 0; i < coefficients.size(); i++)
//      coefficients[i]->interpolate(ufc.w[i], ufc_a.cell, *ufc_a.coefficient_elements[i], *cell);
    
    // Tabulate dofs for each dimension
    for (dolfin::uint i = 0; i < ufc_a.form.rank(); i++)
      ufc_a.dof_maps[i]->tabulate_dofs(ufc_a.dofs[i], ufc_a.mesh, ufc_a.cell);

    for(unsigned int i = 0; i < ufc_a.local_dimensions[0]; i++)
      dofs[i] = map[ (ufc_a.dofs[0])[i] ];

    // Tabulate cell tensor
    ufc_a.cell_integrals[0]->tabulate_tensor(ufc_a.A, ufc_a.w, ufc_a.cell);

    // Add entries to global tensor
    MatSetValues(A, ufc_a.local_dimensions[0], dofs, ufc_a.local_dimensions[0], dofs, ufc_a.A, ADD_VALUES);
//    A.add(ufc_a.A, ufc.local_dimensions, dofs);
    counter++;
  }
  cout << "Finished assembly " << this_process << "  " << toc() << endl;

  cout << "Starting finalise assmebly " << this_process << endl;
  tic();
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  real time2 = toc();
  cout << "Finished finalise assembly " << this_process << "  " << time2 << endl;

  tic();
  Matrix B;
  Assembler assembler;
  assembler.assemble(B, a, mesh);
  real time1 = toc();
  cout << "Standard assemble time " << time1 << endl;

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

  real bc_values[3] = {0, 0, 0};
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
