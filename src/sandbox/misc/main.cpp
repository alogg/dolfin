// Place for random tests

#include <dolfin.h>

using namespace dolfin;

int main()
{
  UnitSquare mesh(10, 10);
  mesh.disp();

  StiffnessMatrix A(mesh);
  A.disp();

  return 0;
}
