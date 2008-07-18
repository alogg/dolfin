#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main(int args, char* argv[])
{
  if(!argv[1])
    error("This program requires an integer command line arugment, e.g. ./demo 100");   

  // Create mesh and forms
  const int N = atoi(argv[1]);

  // Create mesh and forms
  UnitSquare mesh(N, N);
  Function f(mesh,2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  //dolfin::uint NN = (dolfin::uint)sqr(N+1);
  //cout << "Matrix dimension: " << NN << endl; 

  // Assemble uBLAS matrix and vector
  cout << "uBLAS ------------------------------------------------------" << endl;
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

  // Assemble uBLAS matrix and vector
  printf("MTL4 with specified number of nonzeros----------------------\n");
  Assembler ass_mtl4(mesh); 
  MTL4Matrix A_mtl4(A_ublas.size(0), A_ublas.size(1), 5);

  tic(); 
  ass_mtl4.assemble(A_mtl4, a, false); 
  real mtl4_assembly_spec = toc();
  cout << "Assembly time: " << mtl4_assembly_spec << endl;

  tic(); 
  ass_mtl4.assemble(A_mtl4, a, false); 
  real mtl4_reassembly_spec = toc();
  cout << "Re-assembly time: " << mtl4_reassembly_spec << endl;

  // Assemble MTL4 matrix and vector
  printf("MTL4 without specified numer of nonzeros---------------------\n");
  Assembler ass_mtl4b(mesh); 
  MTL4Matrix A_mtl4b(A_ublas.size(0), A_ublas.size(1));

  tic(); 
  ass_mtl4b.assemble(A_mtl4b, a, false); 
  real mtl4_assembly = toc();
  cout << "Assembly time: " << mtl4_assembly << endl;

  tic(); 
  ass_mtl4b.assemble(A_mtl4b, a, false); 
  real mtl4_reassembly = toc();
  cout << "Re-assembly time: " << mtl4_reassembly << endl;

  MTL4Vector b_mtl4b(A_ublas.size(0));
  tic();
  ass_mtl4b.assemble(b_mtl4b, L, false);
  real mtl4_vector_assembly = toc();
  cout << "Vector assembly time: " << mtl4_vector_assembly << endl;


  cout << "-------------------------------------" << endl;
  cout << "Summary:          MTL4 (a), MTL4 (b), uBLAS" << endl;
  cout << "First matrix assembly:   " << mtl4_assembly        << "  " << mtl4_assembly_spec   << "  " << ublas_assembly        << endl;
  cout << "Matrix re-assembly:      " << mtl4_reassembly      << "  " << mtl4_reassembly_spec << "  " << ublas_reassembly      << endl;
  cout << "Vector    assembly:      " << mtl4_vector_assembly << "  " << mtl4_vector_assembly << "  " << ublas_vector_assembly << endl;
  cout << "Matrix/vector dimension: " << A_ublas.size(0) << endl;

  return 0;
}
