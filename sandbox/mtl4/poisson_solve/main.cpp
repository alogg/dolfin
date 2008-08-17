#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main()
{
  // Create mesh and forms
  UnitSquare mesh(128, 128);
  DomainBoundary boundary;
  Function g(mesh, 0.0);
  DirichletBC bc(g, mesh, boundary); 

  Function f(mesh, 2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  Table table;

  // Linear solver
  KrylovSolver solver(bicgstab, ilu);
  //LUSolver solver;
  
  // Assemble and solve using uBLAS
  cout << "uBLAS ------------------------------------------------------" << endl;
  Assembler ass_ublas(mesh);
  uBlasMatrix<ublas_sparse_matrix> A_ublas;
  uBlasVector b_ublas, x_ublas;

  tic(); 
  ass_ublas.assemble(A_ublas, a); 
  ass_ublas.assemble(b_ublas, L); 
  bc.apply(A_ublas, b_ublas, a);
  solver.solve(A_ublas, x_ublas, b_ublas);
  table("uBLAS", "solve") =  toc();

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
  solver.solve(A_mtl4, x_mtl4, b_mtl4);
  table("MTL4", "solve") =  toc();


  //Function u;
  //u.init(mesh, x_mtl4, a, 1);
  //File file("test.pvd");
  //file << u;
  //x_mtl4.disp();

  table.disp();

  return 0;
}
