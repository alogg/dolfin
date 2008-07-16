#include <dolfin.h>

#include "NSEMomentum3D.h"

#include "../DirectAssembler.h"
#include "../MTL_Matrix.h"
#include "../MTL_Vector.h"
  
using namespace dolfin;

void foo(MTL4_sparse_matrix& A, MTL4_vector& b)
{}

int main(int args, char* argv[])
{
  //dolfin_set("output destination", "silent");

  const int N = atoi(argv[1]);
  const int N_rep = atoi(argv[2]);

  UnitCube mesh(N, N, N);
  //mesh.init(mesh.topology().dim());
  Function f(mesh,1.0);
  NSEMomentum3DBilinearForm a(f,f,f,f,f);

  if(1)
    {
      printf("MTL --------------------------------------------------------\n");
      DirectAssembler ass(mesh); 
      const int NN = 3*(N+1)*(N+1)*(N+1);
      MTL_Matrix A(NN,NN,45);

      ass.assemble(A,a,false);

      MTL4_sparse_matrix& mat = A.mat(); 
      MTL4_vector x1(NN, 3.14);
      MTL4_vector x2(NN, 0.0);

      tic();
      for(int i=0; i<N_rep; i++)
	  x2 = mat*x1;

      printf("%d mat-vec mult, mat: %d-by-%d, elapsed %f seconds\n", 
	     N_rep, NN, NN, toc());
    }

  if(1)
    {
      printf("UBLAS ------------------------------------------------------\n");
      Assembler ass(mesh);
      Matrix A;
      tic(); ass.assemble(A,a); printf("assembly time: %f\n", toc());

      uBlasSparseMatrix& AA = A.down_cast<uBlasSparseMatrix>();
      ublas_sparse_matrix& mat = AA.mat();

      const int NN = A.size(0);

      uBlasVector vec1(NN);
      uBlasVector vec2(NN);

      ublas_vector& x1 = vec1.vec();
      ublas_vector& x2 = vec2.vec();
      
      for(int i=0; i<NN; i++)
	x1[i] = 3.14;

      tic();
      for(int i=0; i<N_rep; i++)
	x2 = prod(mat,x1);

      printf("%d mat-vec mult, mat: %d-by-%d, elapsed %f seconds\n", 
	     N_rep, NN, NN, toc());

    }

  return 0;
}
