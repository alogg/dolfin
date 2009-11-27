// Copyright (C) 2009 Andre Massing 
// Licensed under the GNU GPL Version 2.1.
//
// First added:  2008-10-08
// Last changed: 2009-11-22

#include <set>

#include <dolfin/mesh/dolfin_mesh.h> //Should be included before the CGAL Kernel.
#include <math.h>

#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Tetrahedron_3 Tetrahedron_3;


using namespace dolfin;

//CGAL points
typedef std::vector<Point_3>::iterator P3_iterator;
//dolfin points
typedef std::vector<Point>::iterator P_iterator;

typedef std::vector<Segment_3>::iterator Seg_iterator;
typedef std::vector<Triangle_3>::iterator Tri_iterator;
typedef std::vector<Tetrahedron_3>::iterator Tet_iterator;

int main()
{
  Point e0(0.0,0.0,0.0);
  Point e1(1.0,0.0,0.0);
  Point e2(0.0,1.0,0.0);
  Point e3(0.0,0.0,1.0);

  Point a(0.1, 0.0, 0.0);
  Point b(0.0, 0.1, 0.0);
  Point c(0.0, 0.0, 0.1);
  Point d(0.0, 1.0/3.0, 0.0);

  double height = 1.2;
  Point e(1.0/3.0, height, 0.0);
  Point f(2.0/3.0, height, 0.0);
  Point g(0.5, 0.0, 0.0); 

  //points shifted by (2,0,0)
  Point A(3.0, 0.0, 0.0);
  Point B(2.0, 1.0, 0.0);
  Point C(2.0, 0.0, 1.0);
  Point D(2.0, 0.0, 0.0);

  //points shifted by (-1,0,0)
  Point x(-1.0, 0.0, 0.0);
  Point y(-1.0, 1.0, 0.0);
  Point z(-1.0, 0.0, 1.0);

  std::vector<Point> points;
  points.push_back(e1);
  points.push_back(e2);
  points.push_back(e3);
  points.push_back(a);
  points.push_back(b);
  points.push_back(c);
  points.push_back(d);
  points.push_back(e0);
  points.push_back(A);
  points.push_back(B);
  points.push_back(C);
  points.push_back(D);

  UnitCube cube(3,3,2);
  cout <<"Total number of cells in Cube:" << cube.num_cells() <<endl;

  UnitSphere sphere(3);
  cout <<"Total number of cells in Sphere:" << sphere.num_cells() <<endl;
  
  IntersectionOperator io(cube);
//  IntersectionOperator io(cube,"SimpleCartesian"); //the same as before
//  IntersectionOperator io(cube,"ExactPredicates");

  int counter = 0;

  //Point intersection, point by point
  for (P_iterator i = points.begin(); i != points.end(); ++i)
  {
    ++counter;
    std::set<unsigned int> cells;
    io.all_intersected_entities(*i,cells);
    std::cout <<"Intersection with point :" << counter  <<std::endl;
    for (std::set<unsigned int>::const_iterator ic = cells.begin(); ic != cells.end(); ++ic)
    {
      std::cout <<"Id: " <<*ic <<std::endl;
    }
  }

  //Entire point set at all
  std::set<unsigned int> cells;
  io.all_intersected_entities(points,cells);
  cout <<"Intersection with  pointlist:" <<endl;
  for (std::set<unsigned int>::const_iterator ic = cells.begin(); ic != cells.end(); ++ic)
  {
    std::cout <<"Id: " <<*ic <<std::endl;
  }
  
  //Mesh intersection
  cells.clear();
  io.all_intersected_entities(sphere,cells);
  cout <<"Intersection with  sphere:" <<endl;
  for (std::set<unsigned int>::const_iterator ic = cells.begin(); ic != cells.end(); ++ic)
  {
    std::cout <<"Id: " <<*ic <<std::endl;
  }
  
  return 0;
}


