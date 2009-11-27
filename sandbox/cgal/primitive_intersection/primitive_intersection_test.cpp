#include <dolfin/mesh/CGAL_includes.h>

#include <CGAL/AABB_triangle_primitive.h>

#include <vector>
#include <iostream>
#include <string>

//using dolfin::Point;

//typedef CGAL::Simple_cartesian<double> K;

///Definition of the used geometric primitives types.
//typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
//typedef CGAL::Simple_cartesian<double> Kernel;
//typedef K::Point_3                                       Point_3;
//typedef K::Segment_3                                     Segment_3;
//typedef K::Triangle_3                                    Triangle_3;
//typedef K::Tetrahedron_3                                 Tetrahedron_3;
//typedef K::Iso_cuboid_3				   Iso_cuboid_3;
//typedef K::Ray_3					   Ray_3;
//typedef K::Line_3					   Line_3;

//typedef CGAL::Bbox_3					   Bbox_3;


typedef std::vector<Point_3>				      Point_3s;
//typedef std::vector<Point>				      Points;
typedef std::vector<Segment_3>				      Segments;
typedef std::vector<Triangle_3>                               Triangles;
typedef std::vector<Tetrahedron_3>                            Tetrahedrons;

typedef Point_3s::iterator				    Poi3_Iterator;
//typedef Points::iterator				    Poi_Iterator;
typedef Segments::iterator				    Seg_Iterator;
typedef Triangles::iterator				    Tri_Iterator;
typedef Tetrahedrons::iterator				    Tet_Iterator;

//using namespace dolfin;
//typedef dolfin::MeshPrimitive<TriangleCellPrimitive> CellPrimitive;
typedef CGAL::AABB_triangle_primitive<K,Tri_Iterator> CellPrimitive;
typedef CGAL::AABB_traits<K, CellPrimitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;


using std::cout;
using std::endl;
using std::string;

void print_message(const string & s1, const string & s2, bool result)
{
  if (result)
    cout << s1 <<  " and " << s2 << " intersect" << endl;
  else
    cout << s1 <<  " and " << s2 << " DO NOT intersect" << endl;
}

int main()
{
  //Building points
  //Unit points
  Point_3 e0(0.0,0.0,0.0);
  Point_3 e1(1.0,0.0,0.0);
  Point_3 e2(0.0,1.0,0.0);
  Point_3 e3(0.0,0.0,1.0);
  Point_3 e4(1.0,1.0,1.0);

  
  double frac = 1.0/3.0;
  Point_3 a(frac, frac, 0.0);
  Point_3 b(0.0, frac, 0.0);
  Point_3 c(frac, 2*frac, 0.0);
  Point_3 d(0.0, 2*frac, 0.0);

  Point_3 f(0.0, frac, 0.0);
  Point_3 g(-1.0,0.0,0.0);
  Point_3 h(-1.0,1.0,0.0);

  Point_3 x(-1.0, 0.0, 0.0);
  Point_3 y(-1.0, 1.0, 0.0);

  //points shifted by (2,0,0)
  Point_3 A(3.0, 0.0, 0.0);
  Point_3 B(2.0, 1.0, 0.0);
  Point_3 C(2.0, 0.0, 1.0);
  Point_3 D(2.0, 0.0, 0.0);

//  Point de0(0.0,0.0,0.0);
//  Point de1(1.0,0.0,0.0);
//  Point de2(0.0,1.0,0.0);
//  Point de3(0.0,0.0,1.0);
//  Point de4(1.0,1.0,1.0);
  
  //vector of dolfin points to test intersection incoporation the conversion
  //operator.
  
  Segment_3 seg_1(e0,e1);
  Triangle_3 tri_1(e0,e1,e2);
  Tetrahedron_3 tet_1(e0,e1,e2,e3);
  Tetrahedron_3 tet_2(A,B,C,D);
  Iso_cuboid_3 cub_1(e0,e4);
  Bbox_3 box_1(cub_1.bbox());

  string poistring("Point");
  string segstring("Segment");
  string tristring("Triangle");
  string tetstring("Tetrahedron");
  string cubstring("Cuboid");
  string boxstring("Bbox");
  
  std::cout << "Principial intersection capabilities test: " << std::endl;
  //point primitive intersection test
  print_message(poistring,poistring,CGAL::do_intersect(e0,e0)); //FIXED
  print_message(poistring,segstring,CGAL::do_intersect(e0,seg_1)); //FIXED
  print_message(poistring,tristring,CGAL::do_intersect(e0,tri_1)); 
  print_message(poistring,tetstring,CGAL::do_intersect(e0,tet_1)); //FIXED
  print_message(poistring,cubstring,CGAL::do_intersect(e0,cub_1)); //FIXED
  print_message(poistring,boxstring,CGAL::do_intersect(e0,box_1)); //FIXED

  //segment primitive intersection test
//  print_message(segstring,segstring,CGAL::do_intersect(seg_1,seg_1)); //FAILED
  print_message(segstring,tristring,CGAL::do_intersect(seg_1,tri_1)); 
//  print_message(segstring,tetstring,CGAL::do_intersect(seg_1,tet_1)); //FAILED
  print_message(segstring,boxstring,CGAL::do_intersect(seg_1,box_1)); //FIXED

  //triangles primitive intersection test
  print_message(tristring,tristring,CGAL::do_intersect(tri_1,tri_1));
  print_message(tristring,tetstring,CGAL::do_intersect(tri_1,tet_1));
  print_message(tristring,boxstring,CGAL::do_intersect(tri_1,box_1)); //FAILED

  //tetrahedron primitive intersection test
  print_message(tetstring,tetstring,CGAL::do_intersect(tet_1,tet_1)); //FIXED
  print_message(tetstring,boxstring,CGAL::do_intersect(box_1,tet_1)); //FIXED
  print_message(tetstring,boxstring,CGAL::do_intersect(box_1,tet_1)); //FIXED
  print_message(tetstring,boxstring,CGAL::do_intersect(box_1,tet_2)); //FIXED
  
  std::cout << "Additional intersection test: " << std::endl;
  Triangles triangles;
  triangles.push_back(Triangle_3(e0,a,b));
  triangles.push_back(Triangle_3(b,a,c));
  triangles.push_back(Triangle_3(b,c,d));

  Triangle_3 test_tri(b,x,y);

  int trij = 0;
  for (Tri_Iterator i = triangles.begin(); i != triangles.end(); ++i)
  {
    ++trij;
    std::cout <<"Point " <<std::endl;
    if (CGAL::do_intersect(b,*i))
      std::cout <<"intersects with triangle " << trij <<std::endl;
    else
      std::cout <<"DOES NOT intersect with triangle " << trij <<std::endl;
    std::cout <<"Triangle " <<std::endl;
    if (CGAL::do_intersect(test_tri,*i))
      std::cout <<"intersects with triangle " << trij <<std::endl;
    else
      std::cout <<"DOES NOT intersect with triangle " << trij <<std::endl;
  }

  //Now intersection using CGAL AABB_tree
  Tree tree(triangles.begin(), triangles.end());
  std::cout <<"Intersection detection using AABB tre" << std::endl;
  if(tree.do_intersect(test_tri))
    std::cout << std::endl <<"intersection(s) with triangle  test_tri";
  std::cout << std::endl << tree.number_of_intersected_primitives(test_tri)
	    << " intersection(s) with triangle test_tri " << trij << std::endl;


//  triangles.push_back(Triangle_3(e0,a,b));  //
//  triangles.push_back(Triangle_3(e0,e1,b)); //
//  triangles.push_back(Triangle_3(e0,e1,c)); //
//  triangles.push_back(Triangle_3(x,y,z));   //
//  triangles.push_back(Triangle_3(e0,x,y));  //
//  triangles.push_back(Triangle_3(d,x,y));   //
//  triangles.push_back(Triangle_3(e,f,g));   //
//  triangles.push_back(Triangle_3(A,B,C));   //
//  triangles.push_back(Triangle_3(x,y,d));   //
//  triangles.push_back(Triangle_3(x,y,D));   //

//  std::list<Tetrahedron> tetrahedrons;
//  tetrahedrons.push_back(Tetrahedron(a,b,c,d));
//  tetrahedrons.push_back(Tetrahedron(A,B,C,D));
  return 0;
}
