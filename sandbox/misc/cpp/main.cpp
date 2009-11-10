
#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main (int argc, char* argv[])
{
  UnitSquare mesh(1, 1);

  Poisson::FunctionSpace V(mesh);
  Function v(V);
  mesh.refine();
  mesh.refine();

  GenericVector& x = v.vector();
  for (dolfin::uint i = 0; i < V.dim(); i++)
    x.setitem(i, static_cast<double>(i));

  plot(v);

  mesh.refine();

  plot(v);

  return 0;
}
