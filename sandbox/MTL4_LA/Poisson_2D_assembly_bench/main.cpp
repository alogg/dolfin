#include <dolfin.h>

#include "Poisson.h"

#include "../DirectAssembler.h"
#include "../MTL_Matrix.h"
#include "../MTL_Vector.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  //dolfin_set("output destination", "silent");

  const int N = atoi(argv[1]);

  UnitSquare mesh(N, N);
  //mesh.init(mesh.topology().dim());
  Function f(mesh,2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  if(1)
    {
      printf("MTL --------------------------------------------------------\n");
      DirectAssembler ass(mesh); 
      const int NN = (int)sqr(N+1);

      MTL_Matrix A(NN,NN,5);
      tic(); ass.assemble(A,a,false); printf("assembly time: %f\n", toc());
      tic(); ass.assemble(A,a,false); printf("reassembly time: %f\n", toc());

      MTL_Vector b(NN);
      tic();ass.assemble(b,L,false);printf("vector assembly time: %f\n",toc());
    }

  if(1)
    {
      printf("UBLAS ------------------------------------------------------\n");
      Assembler ass(mesh);
      Matrix A;
      tic(); ass.assemble(A,a); printf("assembly time: %f\n", toc());
      tic(); ass.assemble(A,a,false); printf("reassembly time: %f\n", toc());

      Vector b;
      tic(); ass.assemble(b,L); printf("vector assembly time: %f\n", toc());

      if(N < 5)
	{
	  File fmat("mat.xml");
	  File fvec("vec.xml");
	  fmat << A;
	  fvec << b;
	}
    }

  return 0;
}
