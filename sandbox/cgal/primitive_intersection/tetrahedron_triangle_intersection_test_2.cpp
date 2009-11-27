#include <iostream>
#include <list>

#include <CGAL/AABB_tree.h> // must be inserted before kernel
#include <CGAL/AABB_traits.h>
#include "AABB_tetrahedron_primitive.h"
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Simple_cartesian.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/algorithm.h>
#include <vector>


typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
typedef K::Tetrahedron_3 Tetrahedron;

//types
typedef std::list<Tetrahedron>::iterator Tet_iterator;
typedef std::list<Triangle>::iterator Tri_iterator;
typedef std::list<Point>::iterator P_iterator;

//typedef CGAL::AABB_tetrahedron_primitive<K,Tet_iterator> Primitive;
typedef CGAL::AABB_triangle_primitive<K,Tri_iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Tree;

struct My_Point {
  double x,y,z;
  operator Point() const { return Point(x,y,z); }
};

Point shift_point(const My_Point & p)
{
  return p;
}


int main()
{
    Point a(1.0, 0.0, 0.0);
    Point b(0.0, 1.0, 0.0);
    Point c(0.0, 0.0, 1.0);
    Point d(0.0, 0.0, 0.0);

    //points shifted by (2,0,0)
    Point A(3.0, 0.0, 0.0);
    Point B(2.0, 1.0, 0.0);
    Point C(2.0, 0.0, 1.0);
    Point D(2.0, 0.0, 0.0);

    //points shifted by (-1,0,0)
    Point x(-1.0, 0.0, 0.0);
    Point y(-1.0, 1.0, 0.0);
    Point z(-1.0, 0.0, 1.0);

    std::list<Tetrahedron> tetrahedrons;
    std::list<Triangle> triangles;

    tetrahedrons.push_back(Tetrahedron(a,b,c,d));
    tetrahedrons.push_back(Tetrahedron(A,B,C,D));

    triangles.push_back(Triangle(a,b,c));   //intersects 1st
    triangles.push_back(Triangle(A,B,C));   //intersects 2nd
    triangles.push_back(Triangle(x,y,z));   //intersects none
    triangles.push_back(Triangle(x,y,d));   //intersects 1
    triangles.push_back(Triangle(x,y,D));   //intersects 1st and 2nd

    std::list<Point> points;

    points.push_back(a);
    points.push_back(b);
    points.push_back(c);
    points.push_back(d);
    points.push_back(A);
    points.push_back(B);
    points.push_back(C);
    points.push_back(D);

    //construct AABB tree
//    Tree tree(tetrahedrons.begin(),tetrahedrons.end());
    Tree tree(triangles.begin(),triangles.end());

    //counts intersections for each triangle.
//    for (P_iterator i = points.begin(); i != points.end(); ++i)
//    {
//      ++j;
//      if(tree.do_intersect(*i))
//        std::cout << "intersection(s) with point " << j << std::endl;
//      else
//        std::cout << "no intersection with point" << j << std::endl;

//      std::cout << tree.number_of_intersected_primitives(*i)
//        << " intersection(s) with point number " << j << std::endl;
//    }

    //counts intersections for each triangle.
    int j = 0;
    for (Tri_iterator i = triangles.begin(); i != triangles.end(); ++i)
    {
      ++j;
      if(tree.do_intersect(*i))
        std::cout << "intersection(s) with triangle " << j << std::endl;
      else
        std::cout << "no intersection with triangle" << j << std::endl;

//       //computes intersections with segment query
      std::cout << tree.number_of_intersected_primitives(*i)
        << " intersection(s) with triangle number " << j << std::endl;
    }
    std::cout << tree.number_of_intersected_primitives(triangles.front())
        << " intersections(s) with triangle query" << std::endl;
//    std::cout << A;
    return 0;
}
