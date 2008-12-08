// This file is used for testing parallel assembly
//
// To run this demo, make sure to build DOLFIN with the flag customCxxFlags=-DHAS_PARMETIS=1

#include <dolfin.h>

using namespace dolfin;

int main(int argc, char* argv[])
{
  // Read in mesh from XML file in parallel
  File file("unitsquare.xml.gz");
  LocalMeshData data;
  file >> data;

  return 0;
}
