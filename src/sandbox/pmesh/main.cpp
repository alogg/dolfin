// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Magnus Vikstrøm 
// First added:  2007-10-30
// Last changed: 2008-01-09
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>
#include "pPoisson2D.h"
#include "Poisson2D.h"
#include <iostream>
#include <fstream>
#include "getopts.h"
#include "assemble.h"

using namespace dolfin;

void timer(Mesh& mesh, int num_iterations)
{
  dolfin::cout << "Assembling with sequential assembler." << dolfin::endl;
  set("debug level", -1);
  Poisson2DBilinearForm a;
  Matrix A;
  Assembler assembler(mesh);

  tic();
  for(int i=0; i<num_iterations; ++i)
    assembler.assemble(A, a, true);
  printf("Average assemble time: %.3e\n", toc()/static_cast<real>(num_iterations));
}

void p_timer(Mesh& mesh, MeshFunction<dolfin::uint>& partitions, int num_iterations)
{
  dolfin::cout << "Assembling with parallel assembler." << dolfin::endl;
  set("debug level", -1);
  pPoisson2DBilinearForm a;
  Matrix B;
  pAssembler passembler(mesh, partitions);

  tic();
  for(int i=0; i<num_iterations; ++i)
    passembler.assemble(B, a, true);
  printf("Processor %d: Average assemble time: %.3e\n", dolfin::MPI::processNumber(), toc()/static_cast<real>(num_iterations));
}

int main(int argc, char* argv[])
{
  std::string meshfile, partitionsfile, resultfile;
  std::string assembler = "parallel";
  int cells = 200;
  int num_iterations = 1;
  int num_part = dolfin::MPI::numProcesses();
  bool check = false;
  bool sequential = false;

  Options listOpts;
  int switchInt;
  listOpts.addOption("", "mesh", "Mesh File", true);
  listOpts.addOption("", "partitions", "Mesh partitions file", true);
  listOpts.addOption("", "sequential", "Run with sequential assembler", false);
  listOpts.addOption("", "cells", "Number of cells", true);
  listOpts.addOption("", "num_part", "Number of partitions", true);
  listOpts.addOption("", "resultfile", "File to save results in", true);
  listOpts.addOption("", "num_iterations", "Number of times to assemble", true);
  listOpts.addOption("", "check", "Verify assembly result", false);

  if (listOpts.parse(argc, argv))
    while ((switchInt = listOpts.cycle()) >= 0)
    {
      switch(switchInt)
      {
        case 0:
          meshfile = listOpts.getArgs(switchInt);
          break;
        case 1:
          partitionsfile = listOpts.getArgs(switchInt);
          break;
        case 2:
          sequential = true;
          break;
        case 3:
          cells = atoi(listOpts.getArgs(switchInt).c_str());
          break;
        case 4:
          num_part = atoi(listOpts.getArgs(switchInt).c_str());
          break;
        case 5:
          resultfile = listOpts.getArgs(switchInt);
          break;
        case 6:
          num_iterations = atoi(listOpts.getArgs(switchInt).c_str());
          break;
        case 7:
          check = true;
          break;
        default:
          break;
      }
    }

  Mesh mesh;
  MeshFunction<dolfin::uint>* partitions;
  if(meshfile != "")
  {
    dolfin::cout << "Reading mesh from file: " << meshfile << dolfin::endl;
    mesh = Mesh(meshfile);
  }
  else
  {
    printf("Creating UnitSquare(%d,%d)\n", cells, cells);
    mesh = UnitSquare(cells, cells);
  }
  if(partitionsfile != "")
  {
    dolfin::cout << "Reading partitions from file: " << partitionsfile << dolfin::endl;
    partitions = new MeshFunction<dolfin::uint>(mesh, partitionsfile);
  }
  else
  {
    dolfin::cout << "Partitioning mesh into: " << num_part << " partitions" << dolfin::endl;
    partitions = new MeshFunction<dolfin::uint>(mesh);
    mesh.partition(num_part, *partitions);
  }
  if(sequential)
  {
    dolfin::cout << "Assembling with sequential assembler." << dolfin::endl;
    printf("Number of iteration: %d Number of partitions: %d\n", num_iterations, num_part);
    timer(mesh, num_iterations);
  }
  else
  {
    dolfin::cout << "Assembling with parallel assembler." << dolfin::endl;
    printf("Number of iteration: %d Number of partitions: %d\n", num_iterations, num_part);
    p_timer(mesh, *partitions, num_iterations);
  }

  if(check)
  {
    set("debug level", 1);
    dolfin::cout << "Verifying assembly" << dolfin::endl;
    check_assembly(mesh, *partitions);
  }
  return 0;
}
