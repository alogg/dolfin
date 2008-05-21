#include <dolfin.h>
#include "../assembly_tester.h"
#include "NSEMomentum3D.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  dolfin_set("output destination","silent");

  int N = atoi(argv[1]);
  UnitCube mesh(N, N, N);
  
  Function z(mesh,1.0);
  
  NSEMomentum3DBilinearForm a(z,z,z,z,z);

  assembly_tester(mesh,a,"ICNS 3D",3);

  return 0;
}
