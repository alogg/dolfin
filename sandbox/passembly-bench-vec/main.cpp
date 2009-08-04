// Benchmarking parallel assembly (new DOLFIN version)

#include <dolfin.h>
#include "PoissonP1.h"
#include "PoissonP2.h"

using namespace dolfin;

//namespace Poisson = PoissonP1;
namespace Poisson = PoissonP2;

// Source term
class Source : public Function
{
  void eval(double* values, const double* x) const
  {
    values[0] = sin(x[0]);
    values[1] = sin(x[0]);
  }
};

int main()
{
  // Create mesh
  Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");

  // Create function space
  info(mesh.data());
  Poisson::FunctionSpace V(mesh);
  info(mesh.data());

  // Define variational problem
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f;
  L.f = f;
  VariationalProblem problem(a, L);

  // Avoid direct solver for now, seems to break
  problem.parameters("linear_solver") = "iterative";

  // Compute solution
  Function u;
  problem.solve(u);
  //u.vector().disp();
  
  // Debugging
  double norm = u.vector().norm("l2");
  if (dolfin::MPI::process_number() == 0)
    cout << "Norm of solution vector: " << norm << endl;

  // Plot solution
  //plot(u);

  // Save solution in VTK format
  //File file("solution.pvd");
  //file << u;

  return 0;
}
