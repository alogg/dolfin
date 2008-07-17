// Compile with g++ -O3 -DNDEBUG

#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

using namespace mtl;
using namespace itl;

int main(int argc, char* argv[])
{
  const int size = atoi(argv[1]);
  const int N = size * size; 
  typedef compressed2D<double>  matrix_type;

  // Set up a matrix with 5-point-stencil
  matrix_type A(N, N);
  matrix::laplacian_setup(A, size*size, 1);

  A *= -1; // for positive definite

  // Create an ILU(0) preconditioner
  pc::diagonal<matrix_type> P(A);
  //pc::ilu_0<matrix_type> P(A);

  // Set b such that x == 1 is solution; start with x == 0
  dense_vector<double> x(N, 1.0), b(N);
  b = A*x; 
  x= 0;
    
  // Termination criterion: r < 1e-6 * b or N iterations
  noisy_iteration<double> iter(b, 500, 1.e-6);
    
  // Solve Ax == b with preconditioner P
  cg(A, x, b, P, iter);
  bicg(A, x, b, P, iter);

  return 0;
}
