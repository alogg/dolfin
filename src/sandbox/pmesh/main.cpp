// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm 
// First added:  2007-10-30
// Last changed: 2008-01-09
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>
#include <dolfin/pUFC.h>
#include "pPoisson2D.h"
#include "Poisson2D.h"

using namespace dolfin;
//-----------------------------------------------------------------------------
std::string appendRank(std::string base, std::string ext)
{
  std::stringstream stream;
  stream << base << dolfin::MPI::processNumber() << "." << ext;
  return stream.str();
}

void check(Mesh& mesh, MeshFunction<dolfin::uint>& partitions)
{
  // Do normal assembly on process 1
  Poisson2DBilinearForm a;
  pPoisson2DBilinearForm b;
  Matrix A, B;

  if(dolfin::MPI::numProcesses() == 1)
  {
    Assembler assembler(mesh);
    assembler.assemble(A, a, true);
  }
  // Parallel assembly on all procs
  {
    pAssembler passembler(mesh, partitions);
    passembler.assemble(B, b, true);
  }
  pDofMapSet& dof_map_set = b.dofMaps();

  // Would be nice to have automatic testing of B = A * modified dofs
  // Currently just printing so that matrices can be manually inspected
  PetscViewer viewer_a(PETSC_VIEWER_STDOUT_SELF);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "A.m", &viewer_a);
  //PetscViewerSetFormat(viewer_a, PETSC_VIEWER_ASCII_MATLAB);
  MatView(A.mat().mat(), viewer_a);

  dolfin::cout << "Mapping: " << dolfin::endl;
  std::map<const dolfin::uint, dolfin::uint> map = dof_map_set[0].getMap();
  for(dolfin::uint i=0; i<map.size(); ++i)
  {
    dolfin::cout << i << " => " << map[i] << dolfin::endl;
  }
  PetscViewer viewer(PETSC_VIEWER_STDOUT_WORLD);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "B.m", &viewer);
  //PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  MatView(B.mat().mat(), viewer);
}

void timer(Mesh& mesh, MeshFunction<dolfin::uint>& partitions, int num_iterations)
{
  set("debug level", -1);
  Poisson2DBilinearForm a;
  Matrix A;
  Assembler assembler(mesh);
  tic();
  for(int i=0; i<num_iterations; ++i)
    assembler.assemble(A, a, true);
  std::cout << "Average assemble time: " << toc()/num_iterations << std::endl;
}

void p_timer(Mesh& mesh, MeshFunction<dolfin::uint>& partitions, int num_iterations)
{
  set("debug level", -1);
  pPoisson2DBilinearForm a;
  Matrix B;
  pAssembler passembler(mesh, partitions);
  tic();
  for(int i=0; i<num_iterations; ++i)
    passembler.assemble(B, a, true);
  std::cout << "Processor " << dolfin::MPI::processNumber() << " Average assemble time: " << toc()/num_iterations << std::endl;
}


int main(int argc, char* argv[])
{
  if(argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << 
      " <mesh_file> <num_iterations>\n";
    exit(1);
  }

  dolfin_init(argc, argv);

  char* infilename = argv[1];
  int num_iterations = atoi(argv[2]);

  Mesh mesh(infilename);

  MeshFunction<dolfin::uint> partitions;
  mesh.partition(dolfin::MPI::numProcesses(), partitions);

  if(dolfin::MPI::numProcesses() == 1)
    timer(mesh, partitions, num_iterations);
  else
    p_timer(mesh, partitions, num_iterations);

  check(mesh, partitions);
}
