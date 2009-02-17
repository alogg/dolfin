#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

int main()
{  
  mtl::compressed2D<double> A(10, 10);
  mtl::dense_vector<double> x(10);
  mtl::dense_vector<double> f(10);

  itl::basic_iteration<double> iter(f, 500, 1.0e-6);

  itl::pc::ilu_0<mtl::compressed2D<double> > P(A);
  //itl::pc::ic_0<mtl::compressed2D<double> > P(A);

  itl::cg(A, x, f, P, iter);
}













