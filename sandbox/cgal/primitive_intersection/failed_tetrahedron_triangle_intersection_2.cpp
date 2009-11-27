#include <CGAL/Bbox_3.h>
//#include <CGAL/Simple_cartesian.h> 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <iostream>

typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Tetrahedron_3 Tetrahedron_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;
typedef K::Ray_3 Ray_3;
typedef K::Line_3 Line_3;

typedef CGAL::Bbox_3 Bbox_3;

using std::cout;
using std::endl;

int main()
{
  Point_3 A(0.3,0.0,0.0);
  Point_3 B(0.3,0.0,0.1);
  Point_3 C(0.3,0.1,0.1);
  Point_3 D(0.4,0.1,0.1);

  Point_3 a(0.2,0.0,0.0);
  Point_3 b(0.2,0.1,0.0);
  Point_3 c(0.3,0.1,0.1);

  Triangle_3 tri_1(a,b,c);
  Tetrahedron_3 tet_2(A,B,C,D);

  if (CGAL::do_intersect(tet_2,tri_1))
      cout <<"Intersection of Tetraeder  with face triangle" << endl;
  else 
      cout <<"NO Intersection." << endl;

  return 0;
}
