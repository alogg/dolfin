#include <dolfin.h>
#include <math.h>

#include <dolfin/mesh/Point.h>

#include <CGAL/Bbox_3.h>

#include <CGAL/Simple_cartesian.h> //used kernel
typedef CGAL::Simple_cartesian<double> K;

typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Tetrahedron_3 Tetrahedron_3;
typedef K::Iso_cuboid_3 Iso_cuboid_3;
typedef K::Ray_3 Ray_3;
typedef K::Line_3 Line_3;

typedef CGAL::Bbox_3 Bbox_3;

using namespace dolfin;

void print_point(const  Point & p)
{
  std::cout.precision(16);
  std::cout << "<Point x = " << p.x() << " y = " << p.y() << " z = " << p.z() << ">" <<std::endl;
}


int main()
{

  UnitCube cube(3,3,2);
  UnitSphere sphere(3);

  Cell cell_1(cube,102);
  Cell cell_2(sphere, 79);
  
  cout <<"Cube Cell:\n";
  for (VertexIterator v(cell_1); !v.end(); ++v)
  {
    cout <<"vertex  index = " << v->index() <<"\t"
    << "coordinates:\t" << v->point() << endl;
    std::cout <<"\t\t\t\t\t";
    print_point(v->point());
  }

  cout <<"Sphere Cell:\n";
  for (VertexIterator v(cell_2); !v.end(); ++v)
  {
    cout <<"vertex  index = " << v->index() <<"\t"
    << "coordinates:\t" << v->point() << endl;
    std::cout <<"\t\t\t\t\t";
    print_point(v->point());
  }

  if (CGAL::do_intersect(Primitive_Traits<TetrahedronCellPrimitive,K>::datum(cell_1),
			(Primitive_Traits<TetrahedronCellPrimitive,K>::datum(cell_2))))
      cout <<"Cells do intersect." << endl;

  Point_3 a(0.6666666666666666,0.6666666666666666,0.5);
  Point_3 b(1,0.6666666666666666,0.5);
  Point_3 c(1.0,1.0,0.5);
  Point_3 d(1.0,1.0,1.0);

  Point_3 A(0.9037749551350623, 0.9037749551350623, 0.9037749551350623);
  Point_3 B(1.096225044864938, 0.9037749551350623, 0.9037749551350623);
  Point_3 C(1.096225044864938, 0.9037749551350623,1.096225044864938);
  Point_3 D(1.096225044864938, 1.096225044864938, 1.096225044864938);

  Tetrahedron_3 tet_1(a,b,c,d);
  Triangle_3 tri_1(a,b,d);
  Tetrahedron_3 tet_2(A,B,C,D);

  if (CGAL::do_intersect(tet_2, Triangle_3(tet_1[0],tet_1[1],tet_1[2])))
      cout <<"Intersection of Tetraeder 1 with face 1" << endl;
  else 
      cout <<"Intersection of Tetraeder 1 with face 1" << endl;

  if (CGAL::do_intersect(tet_2, Triangle_3(tet_1[0],tet_1[1],tet_1[3])))
      cout <<"Intersection of Tetraeder 1 with face 2" << endl;
  else 
      cout <<"Intersection of Tetraeder 1 with face 2" << endl;

  if (CGAL::do_intersect(tet_2, Triangle_3(tet_1[0],tet_1[2],tet_1[3])))
      cout <<"Intersection of Tetraeder 1 with face 3" << endl;
  else 
      cout <<"Intersection of Tetraeder 1 with face 3" << endl;

  if (CGAL::do_intersect(tet_2, Triangle_3(tet_1[1],tet_1[2],tet_1[3])))
      cout <<"Intersection of Tetraeder 1 with face 4" << endl;
  else 
      cout <<"Intersection of Tetraeder 1 with face 4" << endl;

  return 0;
}
