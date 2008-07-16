#include <dolfin.h>

#include "NSEMomentum3D.h"

#include "../DirectAssembler.h"
#include "../MTL_Matrix.h"
#include "../MTL_Vector.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  //dolfin_set("output destination", "silent");

  const int N = atoi(argv[1]);

  UnitCube mesh(N, N, N);
  //mesh.init(mesh.topology().dim());
  Function f(mesh,1.0);
  NSEMomentum3DBilinearForm a(f,f,f,f,f);
  NSEMomentum3DLinearForm L(f,f,f,f,f,f,f,f);

  if(1)
    {
      printf("MTL --------------------------------------------------------\n");
      DirectAssembler ass(mesh); 
      const int NN = 3*(N+1)*(N+1)*(N+1);
      MTL_Matrix A(NN,NN,45);

      tic(); ass.assemble(A,a,false); printf("assembly time: %f\n", toc());
      tic(); ass.assemble(A,a,false); printf("reassembly time: %f\n", toc());

      MTL_Vector b(NN);
      tic();ass.assemble(b,L,false);printf("vector assembly time: %f\n",toc());
    }

  if(0)
    {
      printf("UBLAS ------------------------------------------------------\n");
      Assembler ass(mesh);
      Matrix A;
      tic(); ass.assemble(A,a); printf("assembly time: %f\n", toc());
      tic(); ass.assemble(A,a,false); printf("reassembly time: %f\n", toc());

      Vector b;
      tic(); ass.assemble(b,L); printf("vector assembly time: %f\n", toc());
    }

  return 0;
}
