#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

int main(int args, char* argv[])
{

  dolfin_set("linear algebra backend","MTL4");

  // Create mesh and forms
  UnitSquare mesh(4, 4);
  Function f(mesh,2.0);
  PoissonBilinearForm a;
  PoissonLinearForm L(f);

  Assembler assembler(mesh);

  //MTL4Matrix A(25,25,5);
  //assembler.assemble(A, a, false); 
  Matrix A;
  assembler.assemble(A, a); 
  
  //A.disp();

  Vector x;
  assembler.assemble(x, L);
  x.disp();

  A *= 2;
  A.disp();

  A /= 2;
  A.disp();

  //MTL4Matrix B=A; // why does this not work?
  Matrix B(A);
  B.disp();

  Vector y;
  A.mult(x, y);
  y.disp();

  double block[] = {10, 20, 30, 40};
  dolfin::uint rows[] = {1, 4};
  dolfin::uint cols[] = {0, 2};

  A.set(block, 2, rows, 2, cols); 
  A.apply();
  A.disp();
  B.disp();

  for(int i=0; i<4; i++) block[i] = 0;

  A.get(block, 2, rows, 2, cols);
  for(int i=0; i<4; i++) cout << block[i] << endl;

  Matrix* C = B.copy();
  B.disp();
  C->set(block, 2, rows, 2, cols);
  C->apply();
  C->disp();

  Array<double> val;
  Array<dolfin::uint> col;
  A.getrow(4,col,val);

  for(dolfin::uint i=0; i<val.size(); i++)
    {
      cout << val[i] << " " << col[i] << endl;
    }

  return 0;
}
