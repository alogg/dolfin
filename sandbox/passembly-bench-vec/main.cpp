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
  Mesh mesh("unitsquare_large.xml.gz");
  //Mesh mesh("unitsquare.xml.gz");
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

  // Assemble matrix
  Matrix A;
  assemble(A, a);

  summary();

  return 0;
}
