
#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main (int argc, char* argv[])
{
  UnitSquare mesh(3, 3);

  Poisson::FunctionSpace V(mesh);
  Function v(V);
  mesh.refine();
  mesh.refine();

  return 0;

  GenericVector& x = v.vector();
  for (dolfin::uint i = 0; i < V.dim(); i++)
    x.setitem(i, static_cast<double>(i));

  plot(v);

  mesh.refine();

  plot(v);

  return 0;
}
