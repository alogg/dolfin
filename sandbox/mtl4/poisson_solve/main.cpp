#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main()
{
  // Create mesh and forms
  UnitSquare mesh(700, 700);
  DomainBoundary boundary;
  Function g(mesh, 0.0);
  DirichletBC bc(g, mesh, boundary); 

  Function f(mesh, 2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  // Assemble and solve using uBLAS
  cout << "uBLAS ------------------------------------------------------" << endl;
  Assembler ass_ublas(mesh);
  uBlasMatrix<ublas_sparse_matrix> A_ublas;
  uBlasVector b_ublas, x_ublas;

  tic(); 
  ass_ublas.assemble(A_ublas, a); 
  ass_ublas.assemble(b_ublas, L); 
  bc.apply(A_ublas, b_ublas, a);
  UmfpackLUSolver lu_ublas;
  lu_ublas.solve(A_ublas, x_ublas, b_ublas);
  real ublas_solve = toc();

  //x_ublas.disp();


  // Assemble and solve using MTL4
  printf("MTL4 with specified number of nonzeros----------------------\n");
  Assembler ass_mtl4(mesh); 
  MTL4Matrix A_mtl4(A_ublas.size(0), A_ublas.size(1), 10);
  MTL4Vector b_mtl4(A_ublas.size(0)), x_mtl4;

  tic();
  ass_mtl4.assemble(A_mtl4, a, false); 
  ass_mtl4.assemble(b_mtl4, L, false);
  bc.apply(A_mtl4, b_mtl4, a);
  //ITLKrylovSolver solver_itl;
  UmfpackLUSolver solver_itl;
  solver_itl.solve(A_mtl4, x_mtl4, b_mtl4);
  real mtl4_solve = toc();


  //Function u;
  //u.init(mesh, x_mtl4, a, 1);
  //File file("test.pvd");
  //file << u;
  //x_mtl4.disp();

  cout << "Solve time (uBLAS): " << ublas_solve << endl;
  cout << "Solve time (MTL4):  " << mtl4_solve << endl;

  return 0;
}
