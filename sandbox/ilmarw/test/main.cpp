#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main()
{ 
  dolfin_set("linear algebra backend", "Assembly");
  
  PoissonBilinearForm a;
  UnitCube mesh(100, 100, 100);
  Matrix A;

  Assembler assembler(mesh);
  assembler.assemble(A, a);

  A.disp();

  return 0;
}
