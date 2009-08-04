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
  Timer t0("MAIN 0: create mesh");
  Mesh mesh("unitsquare_large.xml.gz");
  //Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");
  t0.stop();

  // Create function space

  Timer t1("MAIN 1: create function space");
  Poisson::FunctionSpace V(mesh);
  t1.stop();

  // Define variational problem
  Timer t2("MAIN 2: define variational problem");
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f;
  L.f = f;
  VariationalProblem problem(a, L);
  t2.stop();

  // Avoid direct solver for now, seems to break
  problem.parameters("linear_solver") = "iterative";

  // Assemble matrix
  Timer t3("MAIN 3: assemble");
  Matrix A;
  assemble(A, a);
  t3.stop();

  info("Sum of timings in main: %g", t0.value() + t1.value() + t2.value() + t3.value());
  summary();

  return 0;
}
