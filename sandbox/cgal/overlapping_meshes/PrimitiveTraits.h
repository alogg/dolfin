// =====================================================================================
//
// Copyright (C) 2010-02-02  André Massing
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by André Massing, 2010
//
// First added:  2010-02-02
// Last changed: 2010-02-03
// 
//Author:  André Massing (am), massing@simula.no
//Company:  Simula Research Laboratory, Fornebu, Norway
//
// =====================================================================================

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Point.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Polyhedron_3.h>

namespace dolfin
{

  template <int dim, typename Kernel > struct Primitive_Converter ;

  template<typename Kernel> struct Primitive_Converter<2,Kernel>
  {
    typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
    typedef CGAL::Nef_polyhedron_2<Kernel>  Nef_polyhedron;
    typedef Nef_polyhedron Datum;
    typedef typename Kernel::Point_3 Point_3;

    static const int dim = 2;

//    static Datum datum(const MeshEntity & entity) {

//            return Datum(point);
//    }
  };

  template<typename Kernel> struct Primitive_Converter<3,Kernel>
  {
    typedef CGAL::Polyhedron_3<Kernel>  Polyhedron_3;
    typedef CGAL::Nef_polyhedron_3<Kernel>  Polyhedron;
    typedef Polyhedron Datum;
    typedef typename Kernel::Point_3 Point_3;

    static const int dim = 3;
    static Datum datum(const MeshEntity & entity) {
      VertexIterator v(entity);
      Point_3 p1(v->point());
      ++v;
      Point_3 p2(v->point());
      ++v;
      Point_3 p3(v->point());
      ++v;
      Point_3 p4(v->point());

      Polyhedron_3 Q;
      Q.make_tetrahedron(p1,p2,p3,p4);

      return Datum(Q);
    }
  };

} //end namespace dolfin.
