// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-10-30
// Last changed:
//
// This file is used for testing distribution of the mesh using MPI

#include <dolfin.h>

using namespace dolfin;
//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{

  // Single processor initializes mesh and broadcasts
  UnitCube mesh(9, 9, 9);

  // All processors
  mesh.init();

  int p = MPIManager::processNum();

  std::stringstream meshstream;
  meshstream << "mesh_p" << p << ".xml";
  std::string meshfile = meshstream.str();

  File file(meshfile);
  file << mesh;

  std::stringstream stream;
  stream << "mesh_function_" << p << ".xml";
  std::string filename = stream.str();

  MeshFunction<unsigned int> f(mesh, 0);

  // Single processor partitions mesh and broadcasts
  mesh.partition(6, f);

  // All processors
  File funcfile(filename);
  funcfile << f;

  // This should not be necessary (finalize should be called by ~MPIManager())
  MPIManager::finalize();
  return 0;
}
