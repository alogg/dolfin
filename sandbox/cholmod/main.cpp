
#include <dolfin.h>
#include "Poisson.h"
#include "CholmodCholeskySolver.h"
  
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
  UnitSquare mesh(32, 32);

  dolfin_set("linear algebra backend", "MTL4");

  // Create functions
  Source f(mesh);

  // Create boundary condition
  Function u0(mesh, 0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);
  
  // Define PDE
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  // Linear system
  Matrix A;
  Vector b, x;

  a.updateDofMaps(mesh);
  L.updateDofMaps(mesh);

  assemble_system(A, a.form(), a.coefficients(), a.dofMaps(), 
		  b, L.form(), L.coefficients(), L.dofMaps(),  
		  mesh, bc, 0, 0, 0, true);

  // A.disp();

  // Solve
  CholmodCholeskySolver cs;
  tic();
  cs.solve(A,x,b);
  //solve(A,x,b);
  message("solve time: %f", toc());

  Function u(mesh, x, a);

  // Save solution to file
  //File file("poisson.pvd");
  //file << u;

  return 0;
}
