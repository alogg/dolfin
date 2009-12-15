// Consider moving functionality from here to
// unit tests once it works

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


void test_reconstruction()
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
}

void test_least_squares()
{
  LAPACKMatrix A(4, 3);
  LAPACKVector b(4);

  A(0, 0) = 1;
  A(1, 1) = 1;
  A(2, 2) = 1;
  A(3, 0) = 1;
  A(3, 1) = 1;
  A(3, 2) = 1;

  b[3] = 1;

  info(A, true);
  info(b, true);

  LAPACKSolvers::solve_least_squares(A, b);

  // Should be 0.25, 0.25, 0,25
  info(b, true);
}

int main()
{
  test_least_squares();
  //test_reconstruction();

  return 0;
}
