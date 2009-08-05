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
    //values[1] = sin(x[0]);
  }
};

int main()
{
  // Create mesh
  //Mesh mesh("unitsquare_large.xml.gz");
  Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");

  std::stringstream fname;
  fname << "unitsquare_p" << dolfin::MPI::process_number() << ".xml";
  File outmesh(fname.str());
  outmesh << mesh;

  // Create function space
  Poisson::FunctionSpace V(mesh);

  // Define variational problem
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f;
  L.f = f;
  VariationalProblem problem(a, L);

  // Avoid direct solver for now, seems to break
  problem.parameters("linear_solver") = "iterative";
  problem.parameters["krylov_solver"]("relative_tolerance") = 1.0e-12;

  Function u;
  problem.solve(u);

  double norm = u.vector().norm("l2");
  if (dolfin::MPI::process_number() == 0)
    std::cout << "Norm of solution vector: " << norm << std::endl;

  // Save solution in VTK format
  File file("result/output.pvd");
  file << u;

  summary();

  return 0;
}
