
#include <dolfin.h>
#include "Poisson.h"
  
using namespace dolfin;

int main()
{
  // Source term
  class Source : public Function
  {
  public:
    
    Source(Mesh& mesh) : Function(mesh) {}

    real eval(const real* x) const
    {
      real dx = x[0] - 0.5;
      real dy = x[1] - 0.5;
      return 500.0*exp(-(dx*dx + dy*dy)/0.02);
    }

  };

  // Sub domain for Dirichlet boundary condition
  class DirichletBoundary : public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return x[0] < DOLFIN_EPS && on_boundary;
    }
  };

  // Create mesh
  UnitSquare mesh(512, 512);

  // Create functions
  Source f(mesh);

  // Create boundary condition
  Function u0(mesh, 0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);
  
  // Define PDE
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  // Assembly linear system
  dolfin_set("linear algebra backend", "uBLAS");
  Matrix A;
  Vector b, x;

  assemble(A, a, b, L, bc, mesh);

  Table table("Direct linear solver time");

  // CHOLMOD solver
  CholmodCholeskySolver csolver;
  tic();
  csolver.solve(A, x, b);
  table("CholmodCholeskySolver", "solve time") =  toc();

  // UMFPACK solver
  UmfpackLUSolver usolver;
  tic();
  usolver.solve(A, x, b);
  table("UmfpackLUSolver", "solve time") =  toc();

  table.disp();

  // Save solution to file
  //Function u(mesh, x, a);
  //File file("poisson.pvd");
  //file << u;

  return 0;
}
