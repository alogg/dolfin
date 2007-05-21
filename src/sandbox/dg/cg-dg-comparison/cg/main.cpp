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
//  UnitCube mesh(10, 10, 10);
  UnitCube mesh(5, 5, 5);

  // Create functions
  Source f(mesh);

  FacetNormal n(mesh);
  InvMeshSize h(mesh);

  // Define PDE
//  PoissonBilinearForm a;
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
  File file("poissonCG_stabilised_gamma4_mesh5.pvd");
//  File file("poissonCG_gamma100_mesh10.pvd");
  file << u;

  return 0;
}
