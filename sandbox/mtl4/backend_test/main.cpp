#include <dolfin.h>
#include "Poisson.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  const int N = atoi(argv[1]);
  const int NN = (int)sqr(N+1);
  UnitSquare mesh(N, N);
  Function f(mesh, 1.0);

  PoissonBilinearForm a;
  PoissonLinearForm L(f);
  Assembler ass(mesh);
  
  message("MTL4");
  dolfin_set("linear algebra backend", "MTL4");

  Matrix A(NN,NN);
  ass.assemble(A,a,false);
  
  Vector b(NN);
  ass.assemble(b,L,false);

  return 0;
}
