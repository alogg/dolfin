// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm 
// First added:  2007-10-30
// Last changed: 2008-01-09
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>
#include "Poisson2D.h"

using namespace dolfin;
//-----------------------------------------------------------------------------
std::string appendRank(std::string base, std::string ext)
{
  std::stringstream stream;
  stream << base << dolfin::MPI::processNumber() << "." << ext;
  return stream.str();
}

void check(Mesh& mesh, MeshFunction<dolfin::uint>& partitions, Form& a)
{
  if(dolfin::MPI::numProcesses() == 1)
  {
    Matrix A;
    Assembler assembler(mesh);
    assembler.assemble(A, a, true);

    File file_a("matA.m");
    //file_a << A;
  }
  else
  {
    Matrix B;
    pAssembler passembler(mesh, partitions);
    passembler.assemble(B, a, true);
    B.rename("B", "Parallel matrix B");

    File file_b(appendRank("matB", "m"));

    //file_b << B;
  }
}

void timer(Mesh& mesh, MeshFunction<dolfin::uint>& partitions, Form& a)
{
  Matrix B;
  pAssembler passembler(mesh, partitions);
  tic();
  passembler.assemble(B, a, true);
  std::cout << "Processor " << dolfin::MPI::processNumber() << " assemble time: " << toc() << std::endl;
}

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <num_cells>\n";
    exit(1);
  }

  dolfin_init(argc, argv);

  int num_cells = atoi(argv[1]);

  UnitSquare mesh(num_cells, num_cells);
  MeshFunction<dolfin::uint> partitions;
  mesh.partition(dolfin::MPI::numProcesses(), partitions);

  Poisson2DBilinearForm a;

  //check(mesh, partitions, a);
  timer(mesh, partitions, a);
}
