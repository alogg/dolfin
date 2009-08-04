// Benchmarking parallel assembly (new DOLFIN version)

#include <dolfin.h>
#include "PoissonP1.h"
#include "PoissonP2.h"

using namespace dolfin;

// Source term
class Source : public Function
{
public:
  
  Source(Mesh& mesh) : Function(mesh) {}
  
  void eval(double* values, const double* x) const
  {
    values[0] = sin(x[0]);
  }
  
};

int main()
{
  // Create mesh
  //Mesh mesh("unitsquare_large.xml.gz");
  Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");

  // Define variational problem
  Source f(mesh);
  PoissonP2BilinearForm a;
  PoissonP2LinearForm L(f);
  LinearPDE pde(a, L, mesh);

  // Avoid direct solver for now, seems to break
  dolfin_set("PDE linear solver", "iterative");

  // Solve problem
  Function u;
  pde.solve(u);
  u.vector().disp();

  double norm = u.vector().norm(l2);
  if (dolfin::MPI::processNumber() == 0)
    std::cout << "Norm of solution vector: " << norm << std::endl;

  // Save solution in VTK format
  File file("solution.pvd");
  file << u;

  summary();

  return 0;
}
