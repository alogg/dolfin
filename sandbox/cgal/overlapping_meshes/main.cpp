// =====================================================================================
//
// Copyright (C) 2010-01-28  André Massing
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by André Massing, 2010
//
// First added:  2010-01-28
// Last changed: 2010-01-29
// 
//Author:  André Massing (am), massing@simula.no
//Company:  Simula Research Laboratory, Fornebu, Norway
//
// =====================================================================================

#include <dolfin/mesh/dolfin_mesh.h>
#include <dolfin/common/types.h>
#include <dolfin/plot/dolfin_plot.h>

//#include <dolfin.h>

#include "OverlappingMeshes.h"

using namespace dolfin;
using dolfin::uint;

int main ()
{
//  UnitSquare mesh_1(10,10); 
//  UnitCube mesh_1(10,10,10); 
  UnitCube mesh_1(10,10,10); 

//  UnitSquare mesh_2(10,10,10); 
//  UnitCube mesh_2(10,10,10); 
  UnitSphere mesh_2(10); 

  MeshGeometry& geometry = mesh_2.geometry();
//  double w = 1;
  for (VertexIterator v(mesh_2); !v.end(); ++v)
  {
    double* x = geometry.x(v->index());
//    x[0] = cos(w)*(x[0]-0.5) - sin(w)*(x[1]-0.5);
//    x[1] = sin(w)*(x[0]-0.5) +  cos(w)*(x[1]-0.5);
    x[0] *= 0.5;
    x[1] *= 0.5; 
    x[2] *= 0.5;

    x[0] += 0.4;
    x[1] += 0.4;
    x[2] += 0.4;
  }

  OverlappingMeshes overlap(mesh_1,mesh_2);
  overlap.compute_overlap_map();

  uint_set cells;
  mesh_1.all_intersected_entities(BoundaryMesh(mesh_2),cells);

//  MeshFunction<uint> intersection(mesh_1, mesh_1.topology().dim());
//  intersection = 0;

//  for (uint_set::const_iterator i = cells.begin(); i != cells.end(); i++)
//    intersection[*i] = 1;
//  plot(intersection);

  std::cout <<"overlapped_domain:" << std::endl << overlap.overlapped_domain().str(true);
  plot(overlap.overlapped_domain()); 

  std::cout <<"overlapped_boundary:" << std::endl << overlap.overlapped_boundary().str(true);
  plot(overlap.overlapped_boundary());

  std::cout <<"overlapping_boundary :" << std::endl << overlap.overlapping_boundary().str(true);
  plot(overlap.overlapping_boundary());

  BoundaryMesh boundary_1(mesh_1);
  BoundaryMesh boundary_2(mesh_2);

  cells.clear();
  boundary_1.all_intersected_entities(boundary_2,cells);

  MeshFunction<uint> intersection(boundary_1, boundary_1.topology().dim());
  intersection = 0;

  for (uint_set::const_iterator i = cells.begin(); i != cells.end(); i++)
    intersection[*i] = 1;
  plot(intersection);

  return 0;
}
