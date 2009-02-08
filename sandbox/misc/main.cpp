// Copyright (C) 2009 Marc Spiegelman
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-10-08
// Last changed:  5 Feb 2009 17:16:07
//
// Heavily modified from dolfin/demo/mesh/intersection code
// Program to demonstrated memory leak from  large number of random GTS intersectionDetector tests
// uses boost uniform random number generator for fun
// 
// Note: if you just continually test the point {0.5, 0.5} (comment out the random number generator)
//  the number of intersected cells is 6 rather than 1 and the memory leak is much more pronounced
// 

#include <dolfin.h>
#include <boost/random.hpp>

using namespace dolfin;

typedef boost::mt19937 base_generator_type;

int main()
{
  // Create square mesh
  unsigned int N = 64;
  UnitSquare mesh(N, N);

  //create IntersectionDetector
  IntersectionDetector ID(mesh);
    
  // set initial point
  double x[2] = {0.5, 0.5 };
  Point p(2, x);

  //initialize cell Array
  std::vector<dolfin::uint> cells;
  std::vector<dolfin::uint>::iterator cellIterator;

  //setup boost uniform variate random number generators
  base_generator_type  generator(42u);
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);
  
  dolfin::uint k, k_max = N*N*5000, k_print = 100;
  dolfin::uint n_cells, n_cells_print=2;

  //calculate a random intersection k_max times
  for (k = 0; k < k_max; k++ ) 
  {
    ID.intersection(p, cells);
    n_cells = cells.size();

    //print out information if number of cells in intersection matches n_cells_print (or every k_print grid_sweep)
    if (n_cells >= n_cells_print || k%(N*N*k_print) == 0 ) 
    {
      printf("k=%d, p=[%f, %f], NCells=%d:  ", k, p.x(), p.y(), (int) cells.size()) ;
      for (cellIterator = cells.begin(); cellIterator < cells.end(); cellIterator++) 
      {
        printf("%d ",*cellIterator);
      }
      printf("\n");
    }

    // get new random point (comment out to really blow up)
    p[0] = uni();
    p[1] = uni();

    // clear cells array (cells  grows with every test otherwise, not sure if this generates a memory leak)
    cells.clear(); 
  }
  return(0);
}













