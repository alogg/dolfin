// Testing parallel assembly

#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

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
  //Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");
  //Mesh mesh("unitcube.xml.gz");

  // Uncomment this line to test distribution of built-in meshes
  UnitSquare mesh(2, 2);

  // Create function space
  Poisson::FunctionSpace V(mesh);

  // Define variational problem
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f;
  L.f = f;
  VariationalProblem problem(a, L);

  Function u;
  problem.solve(u);

  return 0;
}
