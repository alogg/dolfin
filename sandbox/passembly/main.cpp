// This file is used for testing parallel assembly

#include <dolfin.h>
#include <dolfin/main/MPI.h>

using namespace dolfin;

int main(int argc, char* argv[])
{
  // Read in mesh from XML file in parallel
  Mesh mesh("unitsquare.xml.gz");
  
  // Store partition to file
  char filename[100];
  sprintf(filename, "mesh_part_%d.xml.gz", dolfin::MPI::process_number());
  File file(filename);
  file << mesh;

  return 0;
}
