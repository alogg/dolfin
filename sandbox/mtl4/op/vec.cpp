#include "boost/tuple/tuple.hpp"
#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main(int args, char* argv[])
{

  dolfin_set("linear algebra backend","MTL4");

  // Create mesh and forms
  UnitSquare mesh(2, 2);
  Function f(mesh,2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  Assembler assembler(mesh);

  MTL4Matrix A(9,9,5);
  assembler.assemble(A, a, false); 

  Vector x;
  assembler.assemble(x, L);
  x.disp();

  //Vector y(x);
  Vector y;
  y = x;
  y.disp();

  x.axpy(2,y);
  x.disp();
  y.disp();

  cout << x.size() << endl;

  real val[9];
  x.get(val);

  for(int i=0; i<9; i++) cout << val[i] << " ";
  cout << endl;

  val[2] = 43;
  x.set(val);
  x.disp();

  x.add(val);
  x.disp();

  real vval[2];
  dolfin::uint cols[2] = {0, 2};
  x.get(vval, 2, cols);

  for(int i=0; i<2; i++) cout << vval[i] << " ";
  cout << endl;

  real dot = x.inner(x);
  cout << dot << endl;

  y = 3.14;
  y.disp();
    
  y *= 2;
  y.disp();

  y /= 2;
  y.disp();

  y += x;
  y.disp();
  
  y -= x;
  y.disp();

  return 0;
}
