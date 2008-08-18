// Copyright (C) 2007-2008 Magnus Vikstrom and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-08-17
// Last changed:

#include <dolfin/common/types.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include "Graph.h"
#include "GraphBuilder.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void GraphBuilder::build(Graph& graph, Mesh& mesh, Graph::Representation rep)
{
  // Clear graph
  graph.clear();

  // Set type
  graph._type = Graph::directed;

  // Build
  if(rep == Graph::dual)
    createMeshDual(graph, mesh);
  else if(rep == Graph::nodal)
    createMeshNodal(graph, mesh);
  else
    error("Graph type unknown");
}
//-----------------------------------------------------------------------------
void GraphBuilder::createMeshNodal(Graph& graph, Mesh& mesh)
{
  error("Partitioning of nodal mesh graphs probably doesn't work. Please test and fix.");

  // Initialise mesh
  mesh.init(0, 0);

  // Get number of vertices, edges and arches
  uint num_vertices = mesh.init(0);
  uint num_edges    = mesh.init(1);
  uint num_arches   = 2*num_edges;

  // Initialise graph
  graph.init(num_vertices, num_edges, num_arches);

  // Create nodal graph. Iterate over edges from all vertices
  uint i = 0, j = 0;
  for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
  {
    graph.vertices[i++] = j;
    uint* entities = vertex->entities(0);
    for (uint k = 0; k < vertex->numEntities(0); k++)
      graph.edges[j++] = entities[k];

    // Replace with this?
    /*
    for (VertexIterator neighbor(vertex); !neighbor.end(); ++neighbor)
    {
      //dolfin_debug1("Edge no %d", j);
		dolfin_debug2("Vertex no %d connected to vertex no %d", vertex->index(), neighbor->index());
      edges[j++] = neighbor->index();
    }
    */
  }
}
//-----------------------------------------------------------------------------
void GraphBuilder::createMeshDual(Graph& graph, Mesh& mesh)
{
  // Initialise mesh
  uint D = mesh.topology().dim();
  mesh.init(D, D);
  MeshConnectivity& connectivity = mesh.topology()(D,D);

  // Get number of vertices, edges and arches
  uint num_vertices = mesh.numCells();
  uint num_arches   = connectivity.size();
  uint num_edges    = num_arches/2;

  // Initialise graph
  graph.init(num_vertices, num_edges, num_arches);

  // Create dual graph. Iterate over neighbors from all cells
  uint i = 0, j = 0;
  for (CellIterator c0(mesh); !c0.end(); ++c0)
  {
    //dolfin_debug1("Cell no %d", i);
    //dolfin_debug1("Cell no %d", c0->index());
    graph.vertices[i++] = j;

    for (CellIterator c1(*c0); !c1.end(); ++c1)
    {
      //dolfin_debug2("Cell no %d connected to cell no %d", c0->index(), c1->index());
      //dolfin_debug2("edges[%d] = %d", j, c1->index());
      graph.edges[j++] = c1->index();

    }
  }
}
//-----------------------------------------------------------------------------
