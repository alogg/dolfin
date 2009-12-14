#include <dolfin.h>
#include "Poisson.h"
#include "P1.h"
#include "P2.h"

using namespace dolfin;

class MyExpression : public Expression
{
  void eval(Array<double>& values, const Array<const double>& x) const
  {
    values[0] = sin(3.0*x[0]);
  }
};

int main()
{
  // Create function spaces
  UnitSquare mesh(4, 4);
  P1::FunctionSpace V(mesh);
  P2::FunctionSpace W(mesh);

  // Create P1 function
  MyExpression f;
  Function v(V);
  v.interpolate(f);

  // Create P2 reconstruction
  Function w(W);
  w.reconstruct(v);

  return 0;
}
