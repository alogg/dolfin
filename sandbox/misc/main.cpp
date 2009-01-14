// Place for random tests

#include <dolfin.h>

using namespace dolfin;

int main()
{
  UnitSquare mesh(2, 2);
  BoundaryMesh boundary(mesh);
  unsigned int a = 2;
  Edge edge(mesh, a);
  cout << "**********OK med edge(mesh, unsigned int)************" << endl;

}
