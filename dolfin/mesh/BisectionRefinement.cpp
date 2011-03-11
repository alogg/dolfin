// Copyright (C) 2006 Johan Hoffman.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2009.
// Modified by Garth N. Wells, 2010-2011.
//
// First added:  2006-11-01
// Last changed: 2011-02-22

#include <dolfin/math/dolfin_math.h>
#include <dolfin/log/dolfin_log.h>
#include "BoundaryMesh.h"
#include "Cell.h"
#include "Edge.h"
#include "Mesh.h"
#include "MeshConnectivity.h"
#include "MeshEditor.h"
#include "MeshFunction.h"
#include "MeshGeometry.h"
#include "MeshTopology.h"
#include "RivaraRefinement.h"
#include "Vertex.h"
#include "BisectionRefinement.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void BisectionRefinement::refine_by_recursive_bisection(Mesh& refined_mesh,
                        const Mesh& mesh, const MeshFunction<bool>& cell_marker)
{
  // Transformation maps
  MeshFunction<dolfin::uint> cell_map;
  std::vector<int> facet_map;

  // Call Rivara refinement
  RivaraRefinement::refine(refined_mesh, mesh, cell_marker, cell_map, facet_map);

  // Transform data
  transform_data(refined_mesh, mesh, cell_map, facet_map);
}
//-----------------------------------------------------------------------------
void BisectionRefinement::transform_data(Mesh& newmesh, const Mesh& oldmesh,
                                         const MeshFunction<uint>& cell_map,
                                         const std::vector<int>& facet_map)
{
  newmesh.data().clear();

  // Rewrite materials
  if (oldmesh.data().mesh_function("material indicators"))
  {
     boost::shared_ptr<MeshFunction<unsigned int> > mat;
    mat = newmesh.data().create_mesh_function("material indicators", newmesh.type().dim());

    for(dolfin::uint i=0; i < newmesh.num_cells(); i++)
      (*mat)[i] = (*oldmesh.data().mesh_function("material indicators"))[cell_map[i]];

    info(TRACE, "MeshData MeshFunction \"material indicators\" transformed.");
  }

  // Rewrite boundary indicators
  if (oldmesh.data().array("boundary facet cells")
        && oldmesh.data().array("boundary facet numbers")
        && oldmesh.data().array("boundary indicators"))
  {

    dolfin::uint num_ent = oldmesh.type().num_entities(0);
    std::vector<dolfin::uint>*  bfc = oldmesh.data().array("boundary facet cells");
    std::vector<dolfin::uint>*  bfn = oldmesh.data().array("boundary facet numbers");
    std::vector<dolfin::uint>*  bi  = oldmesh.data().array("boundary indicators");
    dolfin::uint bi_table_size = oldmesh.num_cells()*num_ent;
    std::vector<int> bi_table;
    bi_table.resize(bi_table_size);

    for(dolfin::uint i= 0 ; i< bi_table_size; i++)
      bi_table[i] = -1;

    for(dolfin::uint i = 0; i < bi->size(); i++)
      bi_table[(*bfc)[i]*num_ent+(*bfn)[i]] = (*bi)[i];

    // Empty loop to count elements
    dolfin::uint bi_size = 0;
    for(dolfin::uint c = 0; c < newmesh.num_cells(); c++)
    {
      for(dolfin::uint f = 0; f < num_ent; f++)
      {
        if (facet_map[c*num_ent+f] != -1)
        {
          dolfin::uint table_map = cell_map[c]*num_ent + facet_map[c*num_ent+f];
          if (bi_table[table_map] != -1)
            bi_size++;
        }
      }
    }

    // Create new MeshData std::vectors for boundary indicators
    std::vector<dolfin::uint>* bfc_new = newmesh.data().create_array("boundary facet cells", bi_size);
    std::vector<dolfin::uint>* bfn_new = newmesh.data().create_array("boundary facet numbers", bi_size);
    std::vector<dolfin::uint>* bi_new  = newmesh.data().create_array("boundary indicators", bi_size);

    // Main transformation loop
    dolfin::uint number_bi = 0;
    for(dolfin::uint c = 0; c < newmesh.num_cells(); c++)
    {
      for(dolfin::uint f = 0; f < num_ent; f++)
      {
        if (facet_map[c*num_ent+f] != -1)
        {
          dolfin::uint table_map = cell_map[c]*num_ent + facet_map[c*num_ent+f];
          if (bi_table[table_map] != -1)
          {
            (*bfc_new)[number_bi] = c;
            (*bfn_new)[number_bi] = f;
            (*bi_new)[number_bi] = bi_table[table_map];
            number_bi++;
          }
        }
      }
    }
    info(TRACE, "MeshData \"boundary indicators\" transformed.");
  }
}
//-----------------------------------------------------------------------------
