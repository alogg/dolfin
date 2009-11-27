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


int main()
{
  unsigned int N_max = 150;
  for (unsigned int N =10; N <= N_max; )
  {
//    UnitCube mesh(N,N,N);
    UnitSquare mesh(N,N);
    IntersectionDetector  old_ic(mesh);
    IntersectionOperator  new_ic(mesh);

    int max = 3000000;

    double X[3]  = {0.0, 0.0, 0.0};

    srand(1);
    
    timer t0;

    std::cout <<"#Starting  intersection with random points using CGAL Interface\t" ;
    unordered_set s_cells;
    for (int i = 1; i <= max; ++i)
    {

      X[0] = std::rand()/static_cast<double>(RAND_MAX);
      X[1] = std::rand()/static_cast<double>(RAND_MAX);
//      X[2] = std::rand()/static_cast<double>(RAND_MAX);

      Point p(X[0],X[1],X[2]);

      new_ic.all_intersected_entities(p,s_cells);
    }

    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;
//    std::cout <<"Lenght of cell container: " <<s_cells.size() <<std::endl;

    srand(1);
    t0.restart();

    std::cout <<"#Starting  intersection with random points using GTS Interface\t";
    std::vector<unsigned int> v_cells;
    for (int i = 1; i <= max; ++i)
    {

      X[0] = std::rand()/static_cast<double>(RAND_MAX);
      X[1] = std::rand()/static_cast<double>(RAND_MAX);
//      X[2] = std::rand()/static_cast<double>(RAND_MAX);

      Point p(X[0],X[1],X[2]);

      old_ic.intersection(p,v_cells);
    }

    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;
//    std::cout <<"Lenght of cell container: " <<v_cells.size() <<std::endl;
    
    N += 10;
  }

  return 0;
}
