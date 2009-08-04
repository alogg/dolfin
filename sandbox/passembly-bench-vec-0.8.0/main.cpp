// Benchmarking parallel assembly (old DOLFIN version / Parafin)

#include <dolfin.h>
#include "Poisson.h"
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
    values[1] = sin(x[0]);
  }
  
  unsigned int rank() const
  {
    return 1;
  }
  
  unsigned int dim(unsigned int i) const
  {
    return 2;
  }
  
};

int main()
{
  // Create mesh
  Mesh mesh("unitsquare_large.xml.gz");
  //Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");

  // Define variational problem
  Source f(mesh);
  // PoissonBilinearForm a;
  // PoissonLinearForm L(f);
  PoissonP2BilinearForm a;
  PoissonP2LinearForm L(f);
  LinearPDE pde(a, L, mesh);

  // Avoid direct solver for now, seems to break
  dolfin_set("PDE linear solver", "iterative");

  // Assemble matrix
  Matrix A;
  assemble(A, a, mesh);

  summary();

  return 0;
}
