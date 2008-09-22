#include "boost/tuple/tuple.hpp"
#include <dolfin.h>
#include "Poisson.h"


using namespace dolfin;


int main(int args, char* argv[])
{
  if(!argv[1])
    error("This program requires an integer command line arugment, e.g. ./demo 100");   

  Table table;

  // Create mesh and forms
  const int N = atoi(argv[1]);

  // Create mesh and forms
  UnitSquare mesh(N, N);
  Function f(mesh,2.0);
  PoissonBilinearForm a;

  Assembler assembler(mesh);

  cout << "Assembling uBLAS" << endl; 
  {
    uBLASSparseMatrix A;
    tic(); 
    assembler.assemble(A, a); 
    table("uBLAS", "assembly") = toc();
  }  

  cout << "Assembling MTL4" << endl; 
  {
    tic(); 
    MTL4Matrix A;
    assembler.assemble(A, a); 
    table("MTL4", "assembly") = toc();
  }  

  cout << "Assembling PETSc" << endl; 
  {
    tic(); 
    PETScMatrix A;
    assembler.assemble(A, a); 
    table("PETSc", "assembly") = toc();
  }  

  cout << "Assembling Trilinos" << endl; 
  {
    tic(); 
    PETScMatrix A;
    assembler.assemble(A, a); 
    table("Trilinos", "assembly") = toc();
  }  

  table.disp();

  return 0;
}
