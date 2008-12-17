// This file is used for testing parallel assembly
//
// To run this demo, make sure to build DOLFIN with the flag customCxxFlags=-DHAS_PARMETIS=1

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

  dolfin::uint process_number = dolfin::MPI::process_number();
  char c[100];
  sprintf(c, "mesh_part_%d.xml.gz", process_number);
  File outfile(c);

  outfile << mesh;

  return 0;
}
