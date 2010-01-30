// =====================================================================================
//
// Copyright (C) 2010-01-16  André Massing
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by André Massing, 2010
//
// First added:  2010-01-16
// back changed: 2010-01-27
// 
//Author:  André Massing (am), massing@simula.no
//Company:  Simula Research Laboratory, Fornebu, Norway
//
// =====================================================================================

#include "OverlappingMeshes.h"

#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/NoDeleter.h>

#include <iostream>

using namespace dolfin;
using dolfin::uint;


  OverlappingMeshes::OverlapData::OverlapData(boost::shared_ptr<const Mesh> mesh) 
: mesh(mesh)
{
  intersected_domain = MeshFunction<uint>(*mesh,mesh->topology().dim());
  intersected_domain = 0;
}

OverlappingMeshes::OverlappingMeshes(const Mesh & mesh_1, const Mesh & mesh_2)
  :
  _overlapped_domain(new OverlapData(reference_to_no_delete_pointer(mesh_1)))
  ,_overlapping_domain(new OverlapData(reference_to_no_delete_pointer(mesh_2)))
  ,_overlapped_boundary(new OverlapData(boost::shared_ptr<const Mesh>(new BoundaryMesh(mesh_1))))
  ,_overlapping_boundary(new OverlapData(boost::shared_ptr<const Mesh>(new BoundaryMesh(mesh_2))))
{}

//OverlappingMeshes::OverlappingMeshes(boost::shared_ptr<const Mesh> mesh_1, boost::shared_ptr<const Mesh> mesh_2)
//  :
//     _overlapped_domain(mesh_1)
//    ,_overlapping_domain(mesh_2)
//    _overlapped_boundary(boost::shared_ptr<const Mesh>(new BoundaryMesh(mesh_1)))
//    ,_overlapping_boundary(boost::shared_ptr<const Mesh>(new BoundaryMesh(mesh_2)))
//{}


void OverlappingMeshes::compute_overlap_map()
{

  const Mesh & mesh_1 = *(_overlapped_domain->mesh);
  const Mesh & mesh_2 = *(_overlapping_domain->mesh);  

  const Mesh & boundary_1 = *(_overlapped_boundary->mesh);
  const Mesh & boundary_2 = *(_overlapping_boundary->mesh);  

  EntityEntitiesMap & cell_cell_map = _overlapped_domain->entity_entities_map[&mesh_2];
  EntityEntitiesMap & facet_1_cell_2_map = _overlapped_boundary->entity_entities_map[&mesh_2];
  EntityEntitiesMap & facet_2_cell_1_map = _overlapping_boundary->entity_entities_map[&mesh_1];

  CellIterator cut_cell(mesh_1);
  CellIterator cut_facet(boundary_2);

  //Step 1: 
  //Intersect boundary of mesh_2 with mesh_1.
  //to get the *partially*
  //intersected cells of mesh_1. 
  //This calculates:
  // a) *partially* intersected cells of mesh1
  // b) the exterior facets of mesh_2 which are (part) of the artificial interface.
  // c) *partially* intersected exterior facets of mesh1 and mesh2.

  for (CellIterator cell(boundary_2); !cell.end(); ++cell)
  {
    mesh_1.all_intersected_entities(*cell, facet_2_cell_1_map[cell->index()]);

    if (facet_2_cell_1_map[cell->index()].empty())
      facet_2_cell_1_map.erase(cell->index());
    //If not empty add cell index and intersecting cell index to the map.
    else
    {
      //Iterate of intersected cell of mesh1, find the overlapping cells of
      //mesh 2 and mark cells in mesh1 as partially overlapped.
      for (EntityListIter cell_iter = facet_2_cell_1_map[cell->index()].begin(); cell_iter != facet_2_cell_1_map[cell->index()].end(); ++cell_iter)
      {
        mesh_2.all_intersected_entities(cut_cell[*cell_iter], cell_cell_map[*cell_iter]);
        _overlapped_domain->intersected_domain[*cell_iter] = 1;
      }

      //Compute partially overlapped boundary cells of mesh1 and mesh2.
      //
      //@remark: Clarify whether it is faster to check first if any and then
      //if compute indeces or just try to compute immediately and erase if the cell
      //index container is empty. Rational: We want to avoid a copy operator
      //(linear). A "any intersection" test should have the same complexity as
      //a "compute all intersection" if no intersection occurs. If a
      //intersection occurrs, we want to compute all intersection anyway.
      //1. Version Compute right away and deleting (delete operation is amortized constant) if empty
      //2. Version Check first and compute if intersected.
      //3. Compute and assign if not empty
      //@remark What is faster? Recompute intersecting cells from mesh2 for the
      //exterior facets in mesh1 or map facet index to cell index, and assign
      //their cell set (which might be bigger as the set, which really only
      //intersects the facet).

      EntityList cut_faces;
      boundary_1.all_intersected_entities(*cell,cut_faces);

      if (!cut_faces.empty())
      {
        _overlapping_boundary->intersected_domain[cell->index()] = 1;

        //Compute for each cut exterior facet in mesh1 the cutting cells in
        //mesh2, mark facet as partially overlapped.
        for (EntityListIter face_iter = cut_faces.begin(); face_iter != cut_faces.end(); ++face_iter)
        {
          mesh_2.all_intersected_entities(cut_facet[*face_iter],facet_1_cell_2_map[*face_iter]);
          _overlapped_boundary->intersected_domain[*face_iter] = 1;
        }
      }
      else
      {
        _overlapping_boundary->intersected_domain[cell->index()] = 2;
      }
    }
  }

  //Step 2:
  //Determine all cells of Mesh 1, which are fully overlapped. This is done by
  //going through the cells, check if they are not partially overlapped and
  //therefore  must then be fully overlapped if  any vertex is intersecting
  //mesh2.
  for (CellIterator cell(mesh_1); !cell.end(); ++cell)
  {
    if (_overlapped_domain->intersected_domain[cell->index()] != 1 && mesh_2.any_intersected_entity(VertexIterator(*cell)->point()) != -1)
      _overlapped_domain->intersected_domain[cell->index()] = 2;
  }

  //Step 3:
  //Determine all cells of the boundary of mesh 1, which are fully overlapped.
  //Same method as in Step 2.
  for (CellIterator cell(boundary_1); !cell.end(); ++cell)
  {
    if (_overlapped_boundary->intersected_domain[cell->index()] != 1 
        && mesh_2.any_intersected_entity(VertexIterator(*cell)->point()) != -1)
      _overlapped_boundary->intersected_domain[cell->index()] = 2;
  }
}

const MeshFunction<uint> & OverlappingMeshes::overlapped_domain() const
{
  return _overlapped_domain->intersected_domain;
}

const MeshFunction<uint> & OverlappingMeshes::overlapping_domain() const
{
  return _overlapping_domain->intersected_domain;
}

const MeshFunction<uint> & OverlappingMeshes::overlapped_boundary() const
{
  return _overlapped_boundary->intersected_domain;
}

const MeshFunction<uint> & OverlappingMeshes::overlapping_boundary() const
{
  return _overlapping_boundary->intersected_domain;
}
