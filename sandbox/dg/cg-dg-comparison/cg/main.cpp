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

  // Create mesh
  UnitCube mesh(8, 8, 8);

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

  // Plot solution
  plot(u);

  // Save solution to file
  File file("poisson.pvd");
  file << u;

  return 0;
}
