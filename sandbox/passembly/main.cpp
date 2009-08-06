// Testing parallel assembly

#include <dolfin.h>
#include "Poisson2DP1.h"
#include "Poisson2DP2.h"
#include "Poisson2DP3.h"
#include "Poisson3DP1.h"
#include "Poisson3DP2.h"
#include "Poisson3DP3.h"

using namespace dolfin;

namespace Poisson = Poisson2DP1;
//namespace Poisson = Poisson2DP2;
//namespace Poisson = Poisson2DP3;
//namespace Poisson = Poisson3DP1;
//namespace Poisson = Poisson3DP2;
//namespace Poisson = Poisson3DP3;

// Source term
class Source : public Function
{
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
  //Mesh mesh("unitcube.xml.gz");

  // Store mesh to VTK
  File mesh_file("partitioned_mesh.pvd");
  mesh_file << mesh;

  // Store mesh to XML
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
  u.update();
  file << u;

  summary();

  return 0;
}
