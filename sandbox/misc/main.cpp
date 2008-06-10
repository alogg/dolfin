// Place for random tests

#include <dolfin.h>

using namespace dolfin;

#include <dolfin.h>
#include "MassMatrix2D.h"

using namespace dolfin;

class MySubDomain : public SubDomain
{
  bool inside(const dolfin::real *x, bool on_boundary) const
  {
    return true;
  }
};

int main()
{
  UnitSquare mesh(3, 3);

  PETScMatrix A;
  MassMatrix2DBilinearForm a;
  MySubDomain subDomain;
  assemble(A, a, mesh, subDomain);

  A.disp();
}
