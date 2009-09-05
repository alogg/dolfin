#include <dolfin.h>

using namespace dolfin;

int main (int argc, char* argv[])
{
  UnitSquare mesh(3, 3);
  mesh.init(1, 2);

  dolfin::cout << mesh << dolfin::endl;

  //info(mesh, true);

  return 0;
}
