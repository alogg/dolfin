#include <CGAL/AABB_tree.h> // must be inserted before kernel
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Simple_cartesian.h>

//#include <dolfin/mesh/MeshPrimitives.h>

//using namespace dolfin;
//typedef Triangle_3 Triangle_3;

//typedef Primitive<2> Triangle_Primitive;
//typedef CGAL::AABB_traits<K,Primitive<2> > AABB_PrimitiveTraits;
//typedef CGAL::AABB_tree<AABB_PrimitiveTraits> Tree;
//typedef Point_3 Point;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Tetrahedron_3 Tetrahedron;
typedef K::Triangle_3 Triangle_3;
typedef K::Point_3 Point_3;

typedef std::list<Triangle_3>::iterator Tri_iterator;
typedef CGAL::AABB_triangle_primitive<K,Tri_iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef std::list<Point_3>::iterator P_iterator;

int main()
{
  //Building points
  Point_3 a(1.0, 0.0, 0.0);
  Point_3 b(0.0, 1.0, 0.0);
  Point_3 c(0.0, 0.0, 1.0);
  Point_3 d(0.0, 0.0, 0.0);

  //points shifted by (2,0,0)
  Point_3 A(3.0, 0.0, 0.0);
  Point_3 B(2.0, 1.0, 0.0);
  Point_3 C(2.0, 0.0, 1.0);
  Point_3 D(2.0, 0.0, 0.0);

  //points shifted by (-1,0,0)
  Point_3 x(-1.0, 0.0, 0.0);
  Point_3 y(-1.0, 1.0, 0.0);
  Point_3 z(-1.0, 0.0, 1.0);

  std::list<Point_3> points;
  points.push_back(a);
  points.push_back(b);
  points.push_back(c);
  points.push_back(d);
  points.push_back(A);
  points.push_back(B);
  points.push_back(C);
  points.push_back(D);
  
  std::list<Triangle_3> triangles;
  triangles.push_back(Triangle_3(a,b,c));   //intersects 1st
  triangles.push_back(Triangle_3(A,B,C));   //intersects 2nd
  triangles.push_back(Triangle_3(x,y,z));   //intersects none
  triangles.push_back(Triangle_3(x,y,d));   //intersects 1
  triangles.push_back(Triangle_3(x,y,D));   //intersects 1st and 2nd

//  std::list<Tetrahedron> tetrahedrons;
//  tetrahedrons.push_back(Tetrahedron(a,b,c,d));
//  tetrahedrons.push_back(Tetrahedron(A,B,C,D));

  Tree tree(triangles.begin(),triangles.end());

//  for (P_iterator i = points.begin(); i != points.end(); ++i)
//  {
//    ++j;
//    if(tree.do_intersect(*i))
//      std::cout << "intersection(s) with point " << j << std::endl;
//    else
//      std::cout << "no intersection with point" << j << std::endl;

//    std::cout << tree.number_of_intersected_primitives(*i)
//      << " intersection(s) with triangle number " << j << std::endl;
//  }
    
//    counts intersections for each triangle.
  int j = 0;
  for (Tri_iterator i = triangles.begin(); i != triangles.end(); ++i)
  {
    ++j;
    if(tree.do_intersect(*i))
      std::cout << "intersection(s) with triangle " << j << std::endl;
    else
      std::cout << "no intersection with triangle" << j << std::endl;

    std::cout << tree.number_of_intersected_primitives(*i)
      << " intersection(s) with triangle number " << j << std::endl;
  }
  return 0;
}
