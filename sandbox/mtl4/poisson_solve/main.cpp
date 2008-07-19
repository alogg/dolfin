#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main()
{
  // Create mesh and forms
  UnitSquare mesh(4, 4);
  DomainBoundary boundary;

  Function f(mesh, 2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  // Assemble and solve using uBLAS
  cout << "uBLAS ------------------------------------------------------" << endl;
  Assembler ass_ublas(mesh);
  uBlasMatrix<ublas_sparse_matrix> A_ublas;
  uBlasVector b_ublas, x_ublas;
  Function g(mesh, 0.0);
  DirichletBC bc(g, mesh, boundary); 
  
  ass_ublas.assemble(A_ublas, a); 
  ass_ublas.assemble(b_ublas, L); 
  bc.apply(A_ublas, b_ublas, a);
  UmfpackLUSolver lu_ublas;
  tic(); 
  lu_ublas.solve(A_ublas, x_ublas, b_ublas);
  real ublas_solve = toc();
  cout << "Solve time: " << ublas_solve << endl;

  x_ublas.disp();

/*
  // Assemble and solve using MTL4
  printf("MTL4 with specified number of nonzeros----------------------\n");
  Assembler ass_mtl4(mesh); 
  MTL4Matrix A_mtl4(A_ublas.size(0), A_ublas.size(1), 5);
  MTL4Vector b_mtl4(A_ublas.size(0)), x_mtl4;

  tic();
  ass_mtl4.assemble(A_mtl4, a, false); 
  ass_mtl4.assemble(b_mtl4, L, false);
  bc.apply(A_mtl4, b_mtl4, a);
  ITLKrylovSolver krylov_itl;
  krylov_itl.solve(A_mtl4, x_mtl4, b_mtl4);


  real mtl4_solve = toc();
  cout << "Solve time: " << mtl4_solve << endl;
  for(dolfin::uint i = 0; i < x_mtl4.size(); ++i)
    std::cout<< x_mtl4.vec()[i] << std::endl;
*/

  return 0;
}
