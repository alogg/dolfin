#include <dolfin.h>

#include "NSEMomentum3D.h"

#include "../DirectAssembler.h"
#include "../MTL_Matrix.h"
#include "../MTL_Vector.h"
  
using namespace dolfin;

int main(int args, char* argv[])
{
  //dolfin_set("output destination", "silent");

  if(!argv[1])
    error("This program requires an integer command line arugment, e.g. ./demo 100");   

  const int N = atoi(argv[1]);

  UnitCube mesh(N, N, N);
  //mesh.init(mesh.topology().dim());
  Array<real> f_array(0.0, 0.0, 1.0);
  Function f_vec(mesh, f_array);
  Function f_scalar(mesh, 1.0);
  NSEMomentum3DBilinearForm a(f_vec, f_scalar, f_scalar, f_scalar, f_scalar);
  NSEMomentum3DLinearForm L(f_vec, f_vec, f_vec, f_scalar, f_scalar, f_scalar, f_scalar, f_scalar);

  printf("MTL --------------------------------------------------------\n");
  DirectAssembler ass_mtl4(mesh); 
  const int NN = 3*(N+1)*(N+1)*(N+1);
  MTL_Matrix A_mtl4(NN,NN,45);

  tic(); 
  ass_mtl4.assemble(A_mtl4, a, false); 
  real mtl4_assembly = toc();
  cout << "Assembly time: " << mtl4_assembly << endl;

  tic(); 
  ass_mtl4.assemble(A_mtl4, a, false); 
  real mtl4_reassembly = toc();
  cout << "Re-assembly time: " << mtl4_reassembly << endl;

  MTL_Vector b_mtl4(NN);
  tic();
  ass_mtl4.assemble(b_mtl4, L, false);
  real mtl4_vector_assembly = toc();
  cout << "Vector assembly time: " << mtl4_vector_assembly << endl;


  printf("UBLAS ------------------------------------------------------\n");
  Assembler ass_ublas(mesh);
  Matrix A_ublas;

  tic(); 
  ass_ublas.assemble(A_ublas, a); 
  real ublas_assembly = toc();
  cout << "Assembly time: " << ublas_assembly << endl;

  tic(); 
  ass_ublas.assemble(A_ublas, a, false); 
  real ublas_reassembly = toc();
  cout << "Re-assembly time: " << ublas_reassembly << endl;

  Vector b_ublas;
  tic(); 
  ass_ublas.assemble(b_ublas, L); 
  real ublas_vector_assembly = toc();

  cout << "Vector assembly time: " << ublas_vector_assembly << endl;
  cout << "Summary:                MTL4,   uBLAS" << endl;
  cout << "First matrix assembly " <<  mtl4_assembly        << "  " << ublas_assembly        << endl;
  cout << "Matrix re-assembly    " <<  mtl4_reassembly      << "  " << ublas_reassembly      << endl;
  cout << "Vector    assembly    " <<  mtl4_vector_assembly << "  " << ublas_vector_assembly << endl;

  return 0;
}
