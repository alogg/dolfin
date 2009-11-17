#include <dolfin.h>
#include "M.h"

using namespace dolfin;

int main (int argc, char* argv[])
{

  UnitInterval mesh(10);
  M::FunctionSpace V(mesh);

  Function orig(V);
  orig.vector() = 1.0;

  Function copy = orig;

  Function copy2(orig);

  return 0;
}
