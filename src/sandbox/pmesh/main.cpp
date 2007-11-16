// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-10-30
// Last changed:
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>
#include <mpi.h>

using namespace dolfin;
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  MPIMeshCommunicator comm;

  UnitSquare mesh(3, 3);
  MeshFunction<unsigned int> f(mesh, 0);

  cout << "Starting mesh distribution" << endl;    
  comm.broadcast(mesh);

  cout << "Starting meshfunction distribution" << endl;    
  comm.broadcast(f);

  return 0;
}
