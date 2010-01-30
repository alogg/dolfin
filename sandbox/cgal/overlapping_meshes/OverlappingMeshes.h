// =====================================================================================
//
// Copyright (C) 2010-01-15  André Massing
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by André Massing, 2010
//
// First added:  2010-01-15
// Last changed: 2010-01-28
// 
//Author:  André Massing (am), massing@simula.no
//Company:  Simula Research Laboratory, Fornebu, Norway
//
// =====================================================================================


#ifndef  __OVERLAPPINGMESHES_H
#define  __OVERLAPPINGMESHES_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>

namespace dolfin
{
  class Mesh;

  typedef std::vector<uint> EntityList;
  typedef std::vector<uint>::const_iterator EntityListIter;
  typedef std::map<uint, EntityList> EntityEntitiesMap;
  typedef std::map<const Mesh *,EntityEntitiesMap> MeshCutEntitiesMap;

  ///This class present a collection of overlapping meshes and provides
  ///functionality to compute cell-cell, cell-facet overlaps.
  class OverlappingMeshes {
    
    ///Helper class to store mesh and corresponding data like intersection maps
    ///and meshfunctions.
    struct OverlapData {
      OverlapData(boost::shared_ptr<const Mesh> mesh);
      boost::shared_ptr<const Mesh> mesh;
      MeshFunction<uint> intersected_domain;
      MeshCutEntitiesMap entity_entities_map;
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
    const boost::shared_ptr<OverlapData> _overlapped_domain;
    const boost::shared_ptr<OverlapData> _overlapping_domain;

    const boost::shared_ptr<OverlapData> _overlapped_boundary;
    const boost::shared_ptr<OverlapData> _overlapping_boundary;

  };


} //end namespace dolfin    

#endif   // ----- #ifndef __OVERLAPPINGMESHES_H  -----
