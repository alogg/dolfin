
#include <cstdlib>
#include <iostream>
#include <climits>

//#include <boost/foreach.hpp>
#include <boost/progress.hpp>

#include <dolfin.h>
#include "Projection.h"

using namespace dolfin;

using boost::timer;
using boost::progress_timer;
using boost::progress_display;
//using std::cout;
//using std::endl;



#ifdef HAS_GTS

class F_1 : public Function
{
public:
  void eval(double* values, const double* x) const
  {
    values[0] = x[0];
  }
};

class F_2 : public Function
{
public:
  void eval(double* values, const double* x) const
  {
    values[0] = x[1];
  }
};

class F_3 : public Function
{
public:
  void eval(double* values, const double* x) const
  {
    values[0] = x[2];
  }
};

class F_4 : public Function
{
public:
  void eval(double* values, const double* x) const
  {
    values[0] = sin(3.0*x[0])*sin(3.0*x[1]);
  }
};

class F_5 : public Function
{
public:
  void eval(double* values, const double* x) const
  {
    values[0] = sin(3.0*x[0])*sin(3.0*x[1])*sin(3.0*x[2]);
  }
};

int main()
{
  
//  bool write = false;
  
//  std::cout << "Should  function eval results be written to stdout ? " << std::endl;
//  std::cin >> write;

  //intialize the random number generator
  srand(1);
  
  unsigned int N_max = 100;
  for (unsigned int N =10; N <= N_max; )
  {
    UnitCube mesh(N,N,N);
    
    Projection::FunctionSpace V(mesh);
    Projection::BilinearForm a(V, V);
    Projection::LinearForm L(V);
    
  //  F_1 f_1;
  //  F_2 f_2;
  //  F_3 f_3;
  //  F_4 f_4;

  //  L.f = f_1;
  //  VariationalProblem pde_1(a, L);
  //  Function Pf_1;
  //  pde_1.solve(Pf_1);

  //  L.f = f_2;
  //  VariationalProblem pde_2(a, L);
  //  Function Pf_2;
  //  pde_2.solve(Pf_2);

  //  L.f = f_3;
  //  VariationalProblem pde_3(a, L);
  //  Function Pf_3;
  //  pde_3.solve(Pf_3);

  //  L.f = f_4;
  //  VariationalProblem pde_4(a, L);
  //  Function Pf_4;
  //  pde_4.solve(Pf_4);

    F_5 f_5;
    L.f = f_5;
    VariationalProblem pde_5(a, L);
    Function Pf_5;
    pde_5.solve(Pf_5);

    int max = 100000;

    double X[3]  = {0.0, 0.0, 0.0};
    double value = 0.0;
    
    std::cout <<"#Starting function evaluation of random points" <<std::endl;
    timer t0;

    for (int i = 1; i <= max; ++i)
    {

      X[0] = std::rand()/static_cast<double>(RAND_MAX);
      X[1] = std::rand()/static_cast<double>(RAND_MAX);
      X[2] = std::rand()/static_cast<double>(RAND_MAX);

      //    cout << "X[0] : " << X[0];
      //    cout << "X[1] : " << X[1];
      //    cout << "X[2] : " << X[2];
      //    f_1.eval(&value, X);
      //    f_2.eval(&value, X);
      //    f_3.eval(&value, X);
      //    f_4.eval(&value, X);
      //    f_5.eval(&value, X);
      //    info("f_5(x) = %g", value);

      Pf_5.eval(&value, X);
      //    info("Pf_5(x) = %g", value);

      //increase number of points by 10
      //    i = 10;

    }

    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;
    
    //Increase mesh size
    N += 10;
  }

  return 0;
}


#else

int main()
{
  info("DOLFIN must be compiled with GTS to run this demo.");
  return 0;
}

#endif

