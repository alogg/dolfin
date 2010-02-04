// =====================================================================================
//
// Copyright (C) 2010-01-15  André Massing
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by André Massing, 2010
//
// First added:  2010-01-15
// Last changed: 2010-02-04
// 
//Author:  André Massing (am), massing@simula.no
//Company:  Simula Research Laboratory, Fornebu, Norway
//
// =====================================================================================

#ifndef  __OVERLAPPINGMESHES_H
#define  __OVERLAPPINGMESHES_H

#include <vector>
#include <map>
#include <utility>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Point.h>

#include "PrimitiveTraits.h"

namespace dolfin
{
  class Mesh;

  typedef std::vector<uint> EntityList;
  typedef std::vector<uint>::const_iterator EntityListIter;
  typedef std::map<uint, EntityList> EntityEntitiesMap;
  typedef std::map<uint, EntityList>::const_iterator EntityEntitiesMapIter;
//  typedef std::map<const Mesh *,EntityEntitiesMap> MeshCutEntitiesMap;


  typedef CGAL::Homogeneous<CGAL::Gmpz> Kernel;
  typedef CGAL::Exact_predicates_inexact_constructions_kernel<CGAL::Gmpz> Kernel;
  typedef CGAL::Polyhedron_3<Kernel>  Polyhedron_3;
  typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
  typedef Nef_polyhedron_3 Polyhedron;

  ///Defines a quadrature rule by a pair of points and weights.
  typedef std::pair<std::vector<Point>,std::vector<double> > Quadrature_Rule;  

  class OverlappingMeshes;


  ///@todo:
  ///@todo Need a function which indicates the end of the iteration.
  ///Think about a better design after the first version is running.
  template <int dim> 
  class OverlappedCell
  {
    typedef Primitive_Converter<3,CGAL::Exact_predicates_exact_constructions_kernel> PT;
    typedef PT::Polyhedron Polyhedron;

  public:
    OverlappedCell(const OverlappingMeshes & overlapping_mesh);
    
    /// Step to the next cut cell and computes the relevant quantities aferwards.
    OverlappedCell<dim> & operator++();

    ///Check if we haved reached the end.
    bool end() { return _map_pos == _map_end_pos; }
    
    ///Return a quadrature rule, meaning a list of points and weights.
//    Quadrature_Rule quadrature_rule();

    ///Returns a presentation of the cut cell.
    Polyhedron polyhedron() { return _polyhedron ; }

    ///Returns a polyhedron as the presentation of the union of the overlapping cells.
    Polyhedron overlapping_polyhedron() { return _overlapping_polyhedron; }

    ///Returns a polyhedron as the presentation of the intersected part.
    Polyhedron overlapped_polyhedron();

  private:

    void compute_polyhedrons();

    const OverlappingMeshes &  _overlapping_mesh;

    EntityEntitiesMapIter _map_pos;
    EntityEntitiesMapIter _map_end_pos;

    //Intersected MeshEntity
    MeshEntityIterator _entity_iter;

    //MeshIterator to iterate over the intersecting cells of the other mesh.
    MeshEntityIterator _overlapping_entity_iter;
   
    //Polydreon, which decodes the intersected part of the cell.
    Polyhedron _polyhedron;

    //Polydreon, which decodes the the polyhedron, which is built up by the
    //overlapping cells belonging to the other mesh.
    Polyhedron _overlapping_polyhedron;

  };

  ///This class present a collection of overlapping meshes and provides
  ///functionality to compute cell-cell, cell-facet overlaps.
  ///@todo Improve design, it is very unflexible concerning extension to more than 2 meshes.
  ///e.g. OverlapData should contain a map, which maps the entity_entities_map to the corresponding
  ///mesh, which the entitylist is refering to.
  class OverlappingMeshes {
    
    ///Helper class to store mesh and corresponding data like intersection maps
    ///and meshfunctions.
    struct OverlapData {
      OverlapData(boost::shared_ptr<const Mesh> mesh);
      boost::shared_ptr<const Mesh> mesh;
      MeshFunction<uint> intersected_domain;

      EntityEntitiesMap entity_entities_map;

//      MeshCutEntitiesMap entity_entities_map;
    };

  public:

    ///Constructor takes a list/vector of meshes. The order of meshes defines
    ///also the "overlapp" priority, i.e. if 2 meshes overlap, the one who
    //appears later in the list actually covers the later one.
    OverlappingMeshes(const Mesh & mesh_1, const Mesh & mesh_2);
//    OverlappingMeshes(boost::shared_ptr<const Mesh> mesh_1, boost::shared_ptr<const Mesh> mesh_2);

    ///Computes the overlap mapping. Mesh2  overlaps mesh1. Computes (for
    ///efficient reasons) in addition the boundary overlaps and the artificial
    ///interface.
    void compute_overlap_map();

    //Return meshfunctions. Only first test, think about better design, for
    //example class like OverlapData or suchlike.
    const MeshFunction<uint> & overlapped_domain() const;
    const MeshFunction<uint> & overlapping_domain() const;
    const MeshFunction<uint> & overlapped_boundary() const;
    const MeshFunction<uint> & overlapping_boundary() const;

  private:
    template<int dim>
    friend class OverlappedCell;
    
    const boost::shared_ptr<OverlapData> _overlapped_domain;
    const boost::shared_ptr<OverlapData> _overlapping_domain;

    const boost::shared_ptr<OverlapData> _overlapped_boundary;
    const boost::shared_ptr<OverlapData> _overlapping_boundary;

  };

template <int dim> 
OverlappedCell<dim>::OverlappedCell(const OverlappingMeshes& overlapping_mesh) 
  : _overlapping_mesh(overlapping_mesh)
  
{
  _map_pos = overlapping_mesh._overlapped_domain->entity_entities_map.begin(); 
  _map_end_pos = overlapping_mesh._overlapped_domain->entity_entities_map.end(); 

  ///Oh herregud, please simplify datastructures!!!!
  _entity_iter = MeshEntityIterator(*(overlapping_mesh._overlapped_domain->mesh),
				    overlapping_mesh._overlapped_domain->mesh->topology().dim());
  _entity_iter[_map_pos->first];

  _overlapping_entity_iter = MeshEntityIterator(*(overlapping_mesh._overlapping_domain->mesh),
						overlapping_mesh._overlapping_domain->mesh->topology().dim());
  _overlapping_entity_iter[_map_pos->second.front()];

}

template <int dim>
OverlappedCell<dim> & OverlappedCell<dim>::operator++()
{
  //Update the entity entities mapping.
  ++_map_pos;
  //Bit clumsy, abusing operator[] to set 
  _entity_iter[_map_pos->first];
  _overlapping_entity_iter[_map_pos->second.front()];

  ///@todo Change this. Use lazy initialization.
  compute_polyhedrons();
  return *this;
}

template <int dim>
void OverlappedCell<dim>::compute_polyhedrons()
{
  Polyhedron cell_polyhedron(PT::datum(_entity_iter[_map_pos->first]));
  _overlapping_polyhedron.clear(); 
  for (EntityListIter index = _map_pos->second.begin(); index != _map_pos->second.end() ; ++index)
  {
    _overlapping_polyhedron += PT::datum(_overlapping_entity_iter[*index]);
  }
  cell_polyhedron -= _overlapping_polyhedron;
}

template <int dim>
typename OverlappedCell<dim>::Polyhedron OverlappedCell<dim>::overlapped_polyhedron()
{
  Polyhedron cell_polyhedron(_entity_iter[_map_pos->first]);
  return cell_polyhedron * _overlapping_polyhedron;  
}


} //end namespace dolfin    

#endif   // ----- #ifndef __OVERLAPPINGMESHES_H  -----
