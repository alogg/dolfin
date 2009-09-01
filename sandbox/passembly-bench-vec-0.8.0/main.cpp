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
  Timer t0("MAIN 0: create mesh");
  Mesh mesh("unitsquare_large.xml.gz");
  //Mesh mesh("unitsquare.xml.gz");
  //Mesh mesh("unitsquare_small.xml.gz");
  //Mesh mesh("unitsquare_reallysmall.xml.gz");
  t0.stop();

  Timer t1("MAIN 1: create function space (not used)");
  t1.stop();

  // Define variational problem
  Timer t2("MAIN 2: define variational problem");
  Source f(mesh);
  // PoissonBilinearForm a;
  // PoissonLinearForm L(f);
  PoissonP2BilinearForm a;
  PoissonP2LinearForm L(f);
  LinearPDE pde(a, L, mesh);
  t2.stop();

  // Avoid direct solver for now, seems to break
  dolfin_set("PDE linear solver", "iterative");

  // Assemble matrix
  Timer t3("MAIN 3: assemble");
  Matrix A;
  assemble(A, a, mesh);
  t3.stop();

  message("Sum of timings in main: %g", t0.value() + t1.value() + t2.value() + t3.value());
  summary();

  return 0;
}
