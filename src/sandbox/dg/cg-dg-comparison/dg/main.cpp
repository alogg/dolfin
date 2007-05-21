#include <dolfin.h>
#include "Poisson.h"
#include "Projection.h"

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
//      real dz = x[2] - 0.5;
//      return 500.0*exp(-(dx*dx + dy*dy + dz*dz)/0.02);
      return 500.0*exp(-(dx*dx + dy*dy)/0.02);
    }

  };

  // Create mesh
//  UnitCube mesh(10, 10, 10);
  UnitCube mesh(5, 5, 5);

  // Create functions
  Source f(mesh);
  FacetNormal n(mesh);
  InvMeshSize h(mesh);

  // Define PDE
  PoissonBilinearForm a(n,h);
  PoissonLinearForm L(f);
  LinearPDE pde(a, L, mesh);

  // Solve PDE
  Function u;
  pde.set("PDE linear solver", "direct");
  pde.solve(u);

  // Project solution onto piecewise linears
  Function uu;
  ProjectionBilinearForm aa;
  ProjectionLinearForm LL(u);
  LinearPDE projection(aa, LL, mesh);
  projection.set("PDE linear solver", "direct");
  projection.solve(uu);

  // Plot solution
  plot(u);
  plot(uu);

  // Save solution to file
  File file("poissonDG_stabilised_gamma24_alpha8_mesh5.pvd");
  file << uu;

  return 0;
}
