#include <dolfin.h>
#include "../assembly_tester.h"
#include "Stokes.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  dolfin_set("output destination","silent");

  int N = atoi(argv[1]);
  UnitSquare mesh(N, N);
    
  Function f(mesh,1.0);
  StokesBilinearForm a(f);

  assembly_tester(mesh,a,"Stokes Stabilized (2D)",3);

  return 0;
}
