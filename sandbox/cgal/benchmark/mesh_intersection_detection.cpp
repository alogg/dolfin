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
  unsigned int N_max = 1000;
  for (unsigned int N =1000; N <= N_max; )
  {
//    UnitCube mesh0(N,N,N);
//    UnitCube mesh1(N,N,N);
    UnitSquare mesh0(N,N);
//    UnitSquare mesh1(100,100);
    UnitCircle mesh1(50);

    MeshGeometry & geo1 = mesh1.geometry();
    // Move and scale circle
    for (VertexIterator vertex(mesh1); !vertex.end(); ++vertex)
    {
      double* x = geo1.x(vertex->index());
      x[0] = 0.2*x[0] + 0.5;
      x[1] = 0.2*x[1] + 0.5;
    }

    std::cout <<"#Creating CGAL intersection operator ...";
    timer t0;
    IntersectionOperator  new_ic(mesh0);
    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;

    unordered_set s_cells;
    
    std::cout <<"#Starting  intersection with mesh1 using CGAL Interface\t" ;
    t0.restart();
    new_ic.all_intersected_entities(mesh1,s_cells);

    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;
    std::cout <<"Lenght of cell container: " <<s_cells.size() <<std::endl;

//    for (unordered_set::const_iterator ic = s_cells.begin(); ic != s_cells.end(); ++ic)
//    {
//      std::cout <<"Id: " <<*ic <<std::endl;
//    }


    std::cout <<"#Creating GTS intersection operator ..." ;
    t0.restart();
    IntersectionDetector  old_ic(mesh0);
    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;

//    std::vector<unsigned int> v_cells;
//    std::cout <<"#Starting  intersection with mesh1 using GTS Interface\t" ;

//    t0.restart();
//    old_ic.intersection(mesh1,v_cells);
//    std::cout <<"Time elapsed: " << t0.elapsed() <<std::endl;

//    std::cout <<"Lenght of cell container: " <<v_cells.size() <<std::endl;

//    for (std::vector<unsigned int>::const_iterator ic = v_cells.begin(); ic != v_cells.end(); ++ic)
//    {
//      std::cout <<"Id: " <<*ic <<std::endl;
//    }

    N += 10;
  }

  return 0;
}
