#include <dolfin.h>
#include "../assembly_tester.h"
#include "Elasticity.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  dolfin_set("output destination","silent");

  dolfin_set("linear algebra backend","PETSc");

  int N = atoi(argv[1]);
  UnitCube mesh(N, N, N);
  
  ElasticityBilinearForm a;

  assembly_tester(mesh,a,"Elasticity 3D",3);

  //  dolfin_set("linear algebra backend","Assembly");

  //  assembly_tester(mesh,a,"Elasticity 3D",3);

  return 0;
}
