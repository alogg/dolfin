// This file is used for testing parallel assembly

#include <dolfin.h>
#include <dolfin/main/MPI.h>

using namespace dolfin;

int main(int argc, char* argv[])
{
  // Read in mesh from XML file in parallel
  File file("unitsquare.xml.gz");
  LocalMeshData data;
  file >> data;

  // Partition mesh
  Mesh mesh;
  MeshPartitioning::partition(mesh, data);
  
  // Store partition to file
  char filename[100];
  sprintf(filename, "mesh_part_%d.xml.gz", dolfin::MPI::process_number());
  File outfile(filename);
  outfile << mesh;

  return 0;
}
