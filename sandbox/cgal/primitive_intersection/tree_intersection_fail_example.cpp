#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_tree.h> // *Must* be inserted before kernel!
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Bbox_3.h>

#include <CGAL/Simple_cartesian.h> //used kernel
typedef CGAL::Simple_cartesian<double> K;

///Definition of the used geometric primitives types.
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Tetrahedron_3 Tetrahedron_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;
typedef K::Ray_3 Ray_3;
typedef K::Line_3 Line_3;

typedef CGAL::Bbox_3 Bbox_3;

typedef std::vector<Triangle_3>                               Triangles;
typedef Triangles::iterator				    Tri_Iterator;

typedef CGAL::AABB_triangle_primitive<K,Tri_Iterator> CellPrimitive;
typedef CGAL::AABB_traits<K, CellPrimitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef Tree::Primitive_id Primitive_id;

int main()
{
  
  double frac = 1.0/3.0;
  Point_3 e0(0.0,0.0,0.0);
  Point_3 a(frac, frac, 0.0);
  Point_3 b(0.0, frac, 0.0);
  Point_3 c(frac, 2*frac, 0.0);
  Point_3 d(0.0, 2*frac, 0.0);
  Point_3 x(-1.0, 0.0, 0.0);
  Point_3 y(-1.0, 1.0, 0.0);

  Triangles triangles;
  triangles.push_back(Triangle_3(e0,a,b));
  triangles.push_back(Triangle_3(b,a,c));
  triangles.push_back(Triangle_3(b,c,d));

  Triangle_3 test_tri(b,x,y);
  
  std::cout <<"Intersection detection via CGAL::do_intersect:" << std::endl;
  int trij = 0;
  int counter = 0;
  for (Tri_Iterator i = triangles.begin(); i != triangles.end(); ++i)
  {
    ++trij;
    if (CGAL::do_intersect(test_tri,*i))
    {
      std::cout <<"Test triangle intersects with triangle " << trij <<std::endl;
      ++counter;
    }
    else
      std::cout <<"DOES NOT intersect with triangle " << trij <<std::endl;
  }
  std::cout <<"------------------------------" <<std::endl;
  std::cout <<counter << " intersection(s)  found." <<std::endl;


  //Now intersection using CGAL AABB_tree
  Tree tree(triangles.begin(), triangles.end());
  std::cout << std::endl 
	    <<"Intersection detection using AABB tree" << std::endl;
  std::list<Primitive_id> primitives;
  tree.all_intersected_primitives(test_tri, std::back_inserter(primitives));
  for (std::list<Primitive_id>::const_iterator i = primitives.begin(); i != primitives.end(); ++i)
    std::cout <<"Test triangle intersects with triangle " << (*i - triangles.begin()) << std::endl;
  std::cout <<"------------------------------" <<std::endl;
  std::cout << tree.number_of_intersected_primitives(test_tri)
	    << " intersection(s) found." << std::endl;

  return 0;
}
