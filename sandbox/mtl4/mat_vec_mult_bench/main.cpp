#include <dolfin.h>

#include "NSEMomentum3D.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  //dolfin_set("output destination", "silent");

  dolfin_assert(argv[1]);
  dolfin_assert(argv[2]);

  if(!argv[1] && !!argv[2])
    error("This program requires an two integer command line arugment, e.g. ./demo 100 10");   

  const int N = atoi(argv[1]);
  const int N_rep = atoi(argv[2]);

  UnitCube mesh(N, N, N);
  Array<real> f_array(0.0, 0.0, 1.0);
  Function f_vec(mesh, f_array);
  Function f_scalar(mesh, 1.0);
  NSEMomentum3DBilinearForm a(f_vec, f_scalar, f_scalar, f_scalar, f_scalar);


  cout << "MTL --------------------------------------------------------" << endl;;
  Assembler ass_mtl4(mesh); 
  const int NN = 3*(N+1)*(N+1)*(N+1);

  MTL4Matrix A_mtl4(NN, NN, 45);
  ass_mtl4.assemble(A_mtl4, a, false);

  mtl4_sparse_matrix& mat_mtl4 = A_mtl4.mat(); 
  mtl4_vector x1_mtl4(NN);
  mtl4_vector x2_mtl4(NN);

  for(int i=0; i<NN; i++)
    x1_mtl4[i] = 3.14;

  tic();
  for(int i = 0; i < N_rep; i++)
	  x2_mtl4 = mat_mtl4*x1_mtl4;
  real mtl_time = toc();

  printf("%d mat-vec mult, mat: %d-by-%d, elapsed %f seconds\n", N_rep, NN, NN, mtl_time);

  cout << "UBLAS ------------------------------------------------------" << endl;
  Assembler ass_ublas(mesh);
  Matrix A_ublas;
  
  ass_ublas.assemble(A_ublas, a); 

  uBlasSparseMatrix& AA = A_ublas.down_cast<uBlasSparseMatrix>();
  ublas_sparse_matrix& mat_ublas = AA.mat();

  uBlasVector vec1_ublas(NN);
  uBlasVector vec2_ublas(NN);

  ublas_vector& x1_ublas = vec1_ublas.vec();
  ublas_vector& x2_ublas = vec2_ublas.vec();
      
  for(int i=0; i<NN; i++)
    x1_ublas[i] = 3.14;

  tic();
  for(int i = 0; i < N_rep; i++)
     ublas::axpy_prod(mat_ublas, x1_ublas, x2_ublas, true);
  real ublas_time = toc();

  printf("%d mat-vec mult, mat: %d-by-%d, elapsed %f seconds\n", N_rep, NN, NN, ublas_time);

  cout << "Matrix-vector mult timing (MTL4, uBLAS): " << mtl_time << "  " << ublas_time << endl;

  return 0;
}
