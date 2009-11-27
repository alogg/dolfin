// Copyright (C) 2009 Andre Massing 
// Licensed under the GNU GPL Version 2.1.
//
// First added:  2008-10-08
// Last changed: 2009-11-22
//
// The IntersectionOperator class relies on a CGAL-based search tree structure
// This file illustrates how to use this search structure in combination with
// CGAL primitives.


#include <dolfin/mesh/dolfin_mesh.h> //Should be included before the CGAL Kernel.

#include <CGAL/AABB_tree.h> // *Must* be inserted before kernel!
#include <CGAL/AABB_traits.h>

#include <CGAL/Simple_cartesian.h> 
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Bbox_3.h>

#include <math.h>

using namespace dolfin;

typedef CGAL::Simple_cartesian<double> K;
//Alternative kernel
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;
typedef K::Triangle_3 Triangle_3;
typedef K::Tetrahedron_3 Tetrahedron_3;

typedef TetrahedronCellPrimitive Primitive;
typedef Primitive_Traits<Primitive,K> PT;
typedef MeshPrimitive<PT> CellPrimitive;

typedef CGAL::AABB_traits<K, CellPrimitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

typedef Tree::Object_and_primitive_id Object_and_primitive_id;
typedef Tree::Primitive_id Primitive_id;

typedef std::vector<Point_3>::iterator P_iterator;
typedef std::vector<Segment_3>::iterator Seg_iterator;
typedef std::vector<Triangle_3>::iterator Tri_iterator;
typedef std::vector<Tetrahedron_3>::iterator Tet_iterator;

//testing intersection with dolfin points.

typedef std::vector<Point>				    Points;
typedef Points::iterator				    Poi_Iterator;

void print_point(const  Point & p)
{
  std::cout.precision(16);
  std::cout << "<Point x = " << p.x() << " y = " << p.y() << " z = " << p.z() << ">" <<std::endl;
}


int main()
{
  //Building points
  //Unit points
  Point_3 e0(0.0,0.0,0.0);
  Point_3 e1(1.0,0.0,0.0);
  Point_3 e2(0.0,1.0,0.0);
  Point_3 e3(0.0,0.0,1.0);

  Point_3 a(0.1, 0.0, 0.0);
  Point_3 b(0.0, 0.1, 0.0);
  Point_3 c(0.0, 0.0, 0.1);
  Point_3 d(0.0, 1.0/3.0, 0.0);
  
  double height = 1.2;
  Point_3 e(1.0/3.0, height, 0.0);
  Point_3 f(2.0/3.0, height, 0.0);
  Point_3 g(0.5, 0.0, 0.0); 

  //points shifted by (2,0,0)
  Point_3 A(3.0, 0.0, 0.0);
  Point_3 B(2.0, 1.0, 0.0);
  Point_3 C(2.0, 0.0, 1.0);
  Point_3 D(2.0, 0.0, 0.0);

  //points shifted by (-1,0,0)
  Point_3 x(-1.0, 0.0, 0.0);
  Point_3 y(-1.0, 1.0, 0.0);
  Point_3 z(-1.0, 0.0, 1.0);

  std::vector<Point_3> points;
  points.push_back(a);
  points.push_back(b);
  points.push_back(c);
  points.push_back(d);
  points.push_back(A);
  points.push_back(B);
  points.push_back(C);
  points.push_back(D);
  
  std::vector<Triangle_3> triangles;
  triangles.push_back(Triangle_3(e0,a,b));
  triangles.push_back(Triangle_3(e0,e1,b)); 
  triangles.push_back(Triangle_3(e0,e1,c)); 
  triangles.push_back(Triangle_3(x,y,z));   
  triangles.push_back(Triangle_3(e0,x,y));  
  triangles.push_back(Triangle_3(d,x,y));   
  triangles.push_back(Triangle_3(e,f,g));   
  triangles.push_back(Triangle_3(A,B,C));   
  triangles.push_back(Triangle_3(x,y,d));   
  triangles.push_back(Triangle_3(x,y,D));   

  std::vector<Tetrahedron_3> tetrahedrons;
  tetrahedrons.push_back(Tetrahedron_3(a,b,c,d));
  tetrahedrons.push_back(Tetrahedron_3(A,B,C,D));

  UnitCube mesh(3,3,3);

  for (MeshEntityIterator check_iter(mesh,3); !check_iter.end(); ++check_iter)
  {
    cout <<"Cell Index: " <<check_iter->index() << endl;
    for (VertexIterator v(*check_iter); !v.end(); ++v)
    { 
      cout <<"vertex  index = " << v->index() <<"\t"
           << "coordinates:\t" << v->point() << endl;
      std::cout <<"\t\t\t\t\t";
      print_point(v->point());
    }
  }

  MeshEntityIterator cell_iter(mesh,3);
    
  //Build the search tree.
  Tree tree(cell_iter,cell_iter.end_iterator());

  //Point intersection
//  int pj = 0;
//  for (P_iterator i = points.begin(); i != points.end(); ++i)
//  {
//    ++pj;
//    if(tree.do_intersect(*i))
//      std::cout << std::endl << "intersection(s) with point " << pj;

//    std::cout << std::endl << tree.number_of_intersected_primitives(*i)
//      << " intersection(s) with point number " << pj << std::endl;

    //Compute the cell ids of the corresponding intersection operations.
//    std::vector<Primitive_id> primitives;
//    tree.all_intersected_primitives(*i,std::back_inserter(primitives));
//    std::cout <<"Cell ids for point " << pj << ":" <<std::endl;
//    for (std::vector<Primitive_id>::iterator id_iter = primitives.begin(); 
//        id_iter != primitives.end(); ++id_iter)
//    {
//      std::cout <<"Id : " <<*id_iter <<std::endl;
//    }
//  }
    
  //counts intersections for each triangle.
//  cout << endl <<"Compute all intersection ids for each triangle: "; 
//  int trij = 0;
//  for (Tri_iterator i = triangles.begin(); i != triangles.end(); ++i)
//  {
//    ++trij;
//    if(tree.do_intersect(*i))
//      std::cout << std::endl <<"intersection(s) with triangle " << trij;
//    std::cout << std::endl << tree.number_of_intersected_primitives(*i)
//      << " intersection(s) with triangle number " << trij << std::endl;
    
    //Compute the cell ids of the corresponding intersection operation.
//    std::vector<Primitive_id> primitives;
//    tree.all_intersected_primitives(*i,std::back_inserter(primitives));
//    std::cout <<"Cell ids for triangle " << trij << ":" <<std::endl;
//    for (std::vector<Primitive_id>::iterator id_iter = primitives.begin(); 
//        id_iter != primitives.end(); ++id_iter)
//    {
//      std::cout <<"Id : " <<*id_iter <<std::endl;
//    }
//  }

  return 0;
}


