#include <dolfin.h>

#include "../assembly_tester.h"

#include "PoissonP1.h"
#include "PoissonP2.h"
#include "PoissonP3.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  dolfin_set("output destination","silent");

  int Nx = atoi(argv[1]);
  UnitSquare mesh(Nx, Nx);
  
  PoissonP1BilinearForm a_p1;
  PoissonP2BilinearForm a_p2;
  PoissonP3BilinearForm a_p3;

  assembly_tester(mesh,a_p1,"Poisson Linear",3);
  assembly_tester(mesh,a_p2,"Poisson Quadratic",3);
  assembly_tester(mesh,a_p3,"Poisson Cubic",3);

  return 0;
}
