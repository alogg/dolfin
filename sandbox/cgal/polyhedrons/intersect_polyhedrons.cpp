#include <CGAL/Simple_cartesian.h> 
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>


#include <iostream>

#include <vector>
#include <list>

//typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Exact_predicates_inexact_constructions_kernel<CGAL::Gmpz> K;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>  Polyhedron;
typedef CGAL::Nef_polyhedron_3<K>  Nef_polyhedron;

//typedef CGAL::Homogeneous<CGAL::Gmpz>  K;
//typedef CGAL::Nef_polyhedron_3<K> Nef_polyhedron_3;

typedef K::Tetrahedron_3 Tetrahedron;
typedef K::Triangle_3 Triangle;
typedef K::Segment_3 Segment;
typedef K::Point_3 Point;

typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef std::list<Triangle>::iterator Tri_iterator;

int main()
{
  
  Point a(1.0, 0.0, 0.0);
  Point b(0.0, 1.0, 0.0);
  Point c(0.0, 0.0, 1.0);
  Point d(0.0, 0.0, 0.0);

  Point A(0.1, 0.0, 0.0);
  Point B(0.0, 0.1, 0.0);
  Point C(0.0, 0.0, 2.0);
  Point D(0.0, 0.0, 0.0);

  //Build a point tetrahedron
  Polyhedron P;
  P.make_tetrahedron();

  Polyhedron Q;
  Q.make_tetrahedron(a, b, c, d);
  Nef_polyhedron NQ(Q);

  Polyhedron F;
  F.make_triangle(A,B,C);
  Nef_polyhedron NF(F);
  
  std::cout << "Exact_predicates_exact_constructions_kernel + SNC_indexed_items"
	    << std::endl
	    << "  allows efficient handling of input "
		 "using floating point coordinates"
	    << std::endl;

  Nef_polyhedron N = NF * NQ;
//  Nef_polyhedron N = NF;
  std ::cout <<"Nef_Polyhedron Tetraeder: " <<std::endl;
  std::cout <<	NQ ;
  std ::cout <<"Nef_Polyhedron Triangle: " <<std::endl;
  std::cout <<	NF ;
  std ::cout <<"Nef_Polyhedron Intersection: " <<std::endl;
  std::cout <<	N ;

//  if(N.is_simple()) {
//    N.convert_to_polyhedron(P);
//    std::cout << P;
//  } 
//  else {
//    std::cout << N;
//  }

  
//  CGAL::set_ascii_mode( std::cout);
//  for (Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
//    std::cout << v->point() << std::endl;
//  for (Vertex_iterator v = Q.vertices_begin(); v != Q.vertices_end(); ++v)
//    std::cout << v->point() << std::endl;

//  if (Q.is_closed())

//  Nef_polyhedron_3 NF(F.vertices_begin(),F.vertices_end());

//  Nef_polyhedron_3 N = NF * NQ;
//  std::cout << N;

  return 0;
}
