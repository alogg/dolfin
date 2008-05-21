#include <dolfin.h>
#include "../assembly_tester.h"
#include "Stokes.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  dolfin_set("output destination","silent");

  int N = atoi(argv[1]);
  UnitSquare mesh(N, N);
    
  StokesBilinearForm a;

  assembly_tester(mesh,a,"Stokes TH (2D)",3);

  return 0;
}
