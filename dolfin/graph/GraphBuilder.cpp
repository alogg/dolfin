// Copyright (C) 2007-2008 Magnus Vikstrom and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-08-17
// Last changed:

#include <algorithm>
#include <numeric>
#include <set>
#include <vector>

#include <dolfin/log/log.h>
#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/LocalMeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include "Graph.h"
#include "GraphBuilder.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void GraphBuilder::build(LocalMeshData& mesh_data)
{
  warning("GraphBuilder is highly experimental.");

  std::vector<std::vector<uint> >& cell_vertices = mesh_data.cell_vertices;
  //const std::vector<uint>& global_cell_indices   = mesh_data.global_cell_indices;

  const uint num_gobal_cells    = mesh_data.num_global_cells;
  const uint num_local_cells    = mesh_data.global_cell_indices.size();
  const uint topological_dim    = mesh_data.tdim;
  const uint num_cell_facets    = topological_dim + 1;
  const uint num_facet_vertices = topological_dim;
  const uint num_cell_vertices  = topological_dim + 1;

  // Sort cell vertex indices
  std::vector<std::vector<uint> >::iterator vertices;
  for (vertices = cell_vertices.begin(); vertices != cell_vertices.end(); ++vertices)
    std::sort(vertices->begin(), vertices->end());

  // FIXME: Use std::set or std::vector?
  // Create graph (cell are graph vertices, cell-cell connections are graph edges) 
  std::vector<std::set<uint> > graph(num_local_cells);

  // Get number of cells on each process
  std::vector<uint> cells_per_process = MPI::gather(num_local_cells);

  // Compute offset for going from local to (internal) global numbering
  std::vector<uint> process_offsets(MPI::num_processes());
  for (uint i = 0; i < MPI::num_processes(); ++i)
    process_offsets[i] = std::accumulate(cells_per_process.begin(), cells_per_process.begin() + i, 0); 
  const uint process_offset = process_offsets[MPI::process_number()];

  // Compute local edges (cell-cell connections) using global (internal) numbering
  compute_connectivity(cell_vertices, num_cell_facets, num_facet_vertices, process_offset, graph);

  ///-----------------------------------------------
  /// The rest only applies when running in parallel
  ///-----------------------------------------------

  // Determine candidate ghost cells (graph ghost vertices)
  std::vector<uint> candidate_ghost_cells;
  for (uint i = 0; i < num_local_cells; ++i)
  {
    assert(i < graph.size());
    if (graph[i].size() != num_cell_facets)
      candidate_ghost_cells.push_back(i);
  }  
  cout << "Number of possible boundary cells " << candidate_ghost_cells.size() << endl; 

  // Get number of possible ghost cells coming from each process
  std::vector<uint> ghost_cells_per_process = MPI::gather(candidate_ghost_cells.size());
  
  // Pack local data for candidate ghost cells (global cell index and vertices)
  std::vector<uint> connected_cell_data;
  for (uint i = 0; i < candidate_ghost_cells.size(); ++i)
  {
    // Global (internal) cell index
    connected_cell_data.push_back(candidate_ghost_cells[i] + process_offset);

    // Candidate cell vertices
    std::vector<uint>& vertices = cell_vertices[candidate_ghost_cells[i]];
    for (uint j = 0; j < num_cell_vertices; ++j)
      connected_cell_data.push_back(vertices[j]);  
  }    

  // Prepare package to send (do not send data belonging to this process)
  std::vector<uint> partition;
  std::vector<uint> transmit_data;
  for (uint i = 0; i < MPI::num_processes(); ++i)
  {
    if(i != MPI::process_number())
    {
      transmit_data.insert(transmit_data.end(), connected_cell_data.begin(), 
                           connected_cell_data.end());
      partition.insert(partition.end(), connected_cell_data.size(), i);
    }
  }

  // Set number of candidate ghost cells on this process to zero (not communicated to self)
  ghost_cells_per_process[MPI::process_number()] = 0;

  // Distribute data to all processes
  MPI::distribute(transmit_data, partition);

  // Data structures for unpacking data
  std::vector<std::vector<std::vector<uint> > > candidate_ghost_cell_vertices(MPI::num_processes());
  std::vector<std::vector<uint> > candidate_ghost_cell_global_indices(MPI::num_processes());

  // Unpack data
  uint _offset = 0;
  for (uint i = 0; i < MPI::num_processes()-1; ++i)
  {
    const uint p = partition[_offset]; 
    const uint data_length = (num_cell_vertices+1)*ghost_cells_per_process[p];

    std::vector<uint>& _global_cell_indices         = candidate_ghost_cell_global_indices[p];
    std::vector<std::vector<uint> >& _cell_vertices = candidate_ghost_cell_vertices[p];

    // Loop over data for each cell
    for (uint j = _offset; j < _offset + data_length; j += num_cell_vertices+1)
    {
      assert(partition[j] == p);

      // Get cell global index
      _global_cell_indices.push_back(transmit_data[j]);

      // Get cell vertices
      std::vector<uint> vertices;
      for (uint k = 0; k < num_cell_vertices; ++k)
        vertices.push_back(transmit_data[(j+1)+k]);
      _cell_vertices.push_back(vertices);       
    }
    
    // Update offset
    _offset += data_length;
  }

  // Add off-process (ghost) edges (cell-cell) connections to graph
  std::set<uint> ghost_cell_global_indices;
  for (uint i = 0; i < candidate_ghost_cell_vertices.size(); ++i)
  {
    compute_connectivity(cell_vertices, candidate_ghost_cell_vertices[i], 
                                               candidate_ghost_cell_global_indices[i], 
                                               num_cell_facets, num_facet_vertices, 
                                               graph, ghost_cell_global_indices);
  }

  cout << "Create SCOTCH data" << endl;
  compute_scotch_data(graph, ghost_cell_global_indices, num_gobal_cells);
  cout << "End create SCOTCH data" << endl;
}
//-----------------------------------------------------------------------------
void GraphBuilder::compute_connectivity(const std::vector<std::vector<uint> >& cell_vertices,
                                        uint num_cell_facets, 
                                        uint num_facet_vertices, uint offset,
                                        std::vector<std::set<uint> >& graph)
{
  std::vector<uint>::iterator it;

  // Forward step
  for (uint i = 0; i < cell_vertices.size() - 1; ++i)
  {
    for (uint j =  cell_vertices.size()-1; j > i; --j)
    {
      assert(i != j);

      // Find numer of vertices shared by cells i and j
      std::vector<uint> intersection(num_cell_facets);
      it = std::set_intersection(cell_vertices[i].begin(), cell_vertices[i].end(), 
                                 cell_vertices[j].begin(), cell_vertices[j].end(), 
                                 intersection.begin());
      const uint num_shared_vertices = it - intersection.begin();  

      // Insert edge if cells are neighbours
      if ( num_shared_vertices == num_facet_vertices )
      {
        graph[i].insert(j + offset);
        //cout << "Edge inserted" << endl;
      }
      else if ( num_shared_vertices > num_facet_vertices)
        error("Too many shared vertices. Cannot construct dual graph.");
    }
  }
  // Reverse step (tranpose)
  for (uint i = cell_vertices.size() - 1; i > 0; --i)
  {
    for (uint j = 0; j < i; ++j)
    {
      assert(i != j);

      // Find numer of vertices shared by cell0 and cell1
      std::vector<uint> intersection(num_cell_facets);
      it = std::set_intersection(cell_vertices[i].begin(), cell_vertices[i].end(), 
                                 cell_vertices[j].begin(), cell_vertices[j].end(), 
                                 intersection.begin());
      const uint num_shared_vertices = it - intersection.begin();  

      // Insert edge if cells are neighbours
      if ( num_shared_vertices == num_facet_vertices )
      {
        graph[i].insert(j + offset );
        //cout << "Edge inserted" << endl;
      }
      else if ( num_shared_vertices > num_facet_vertices)
        error("Too many shared vertices. Cannot construct dual graph.");
    }
  }  
}
//-----------------------------------------------------------------------------
dolfin::uint GraphBuilder::compute_connectivity(const std::vector<std::vector<uint> >& cell_vertices,
                                        const std::vector<std::vector<uint> >& candidate_ghost_vertices,
                                        const std::vector<uint>& candidate_ghost_global_indices,
                                        uint num_cell_facets,
                                        uint num_facet_vertices,
                                        std::vector<std::set<uint> >& graph,
                                        std::set<uint>& ghost_cells)
{
  // FIXME: This function can be made more efficient. For example, loop over local
  //        candidate cells (not all cells) and ghost candidates only. 

  // FIXME: It can be made efficient when the SCOTCH numbering is figured out

  const uint num_ghost_vertices_0 = ghost_cells.size();

  std::vector<uint>::iterator it;
  for (uint i = 0; i < cell_vertices.size(); ++i)
  {
    for (uint j = 0; j < candidate_ghost_vertices.size(); ++j)
    {
      // Find numer of vertices shared by cells i and j
      std::vector<uint> intersection(num_cell_facets);
      it = std::set_intersection(cell_vertices[i].begin(), cell_vertices[i].end(), 
                                 candidate_ghost_vertices[j].begin(), candidate_ghost_vertices[j].end(), 
                                 intersection.begin());
      const uint num_shared_vertices = it - intersection.begin();  

      if ( num_shared_vertices == num_facet_vertices )
      {
        graph[i].insert(candidate_ghost_global_indices[j]);
        ghost_cells.insert(candidate_ghost_global_indices[j]);
      }
      else if ( num_shared_vertices > num_facet_vertices)
        error("Too many shared vertices. Cannot construct dual graph.");
    }
  }
  
  // Return number of newly added ghost vertices
  return ghost_cells.size() - num_ghost_vertices_0;
}
//-----------------------------------------------------------------------------
void GraphBuilder::compute_scotch_data(const std::vector<std::set<uint> >& graph,
                               const std::set<uint>& ghost_cells,
                               uint num_global_vertices)
{
  // C-style array indexing
  //const uint baseval = 0;

  // Local data ---------------------------------

  // Number of local graph vertices (cells)
  const uint vertlocnbr = graph.size();
  cout << "Number of local verticies: " << vertlocnbr << endl;

  // Number of local + ghost graph vertices (cells)
  const uint vertgstnbr = graph.size() + ghost_cells.size();
  cout << "Number of edges: " << vertgstnbr << endl;

  // Number of local edges + edges connecting to ghost vertices
  uint edgelocnbr = 0;  
  std::vector<std::set<uint> >::const_iterator vertex;
  for(vertex = graph.begin(); vertex != graph.end(); ++vertex)
    edgelocnbr += vertex->size();
  cout << "Number of local edges: " << edgelocnbr << endl;

  // Local graph layout

  std::vector<uint> vertloctab;
  std::vector<uint> edgeloctab;
  vertloctab.push_back(0);
  for (uint i = 0; i < graph.size(); ++i)
  {
    const std::set<uint>& vertices = graph[i];
    vertloctab.push_back(vertloctab[i] + vertices.size());
    edgeloctab.insert( edgeloctab.end(), vertices.begin(), vertices.end() );
  }

  // Global data ---------------------------------

  // Total  (global) number of vertices (cells) in the graph
  const uint vertglbnbr = num_global_vertices;
  cout << "Number of global cells: " << vertglbnbr << endl;

  // Total (global) number of edges (cell-cell connections) in the graph
  std::vector<uint> num_global_edges = MPI::gather(edgelocnbr);
  const uint edgeglbnbr = std::accumulate(num_global_edges.begin(), num_global_edges.end(), 0);
  cout << "**Number of global edges: " << edgeglbnbr << endl;

  // Number of processes
  const uint procglbnbr = MPI::num_processes();
  cout << "Number of processes: " << procglbnbr << endl;

  // Array containing the number of local vertices (cells) on each process
  std::vector<uint> proccnttab = MPI::gather(graph.size());
  for (uint i = 0; i < MPI::process_number(); ++i)
    cout << "Testing proccnttab " << proccnttab[i] << endl;

  // Array containing . . . . 
  std::vector<uint> procvrttab(MPI::num_processes() + 1);
  for (uint i = 0; i < MPI::num_processes(); ++i)
    procvrttab[i] = std::accumulate(proccnttab.begin(), proccnttab.begin() + i, 0);
  procvrttab[MPI::num_processes()] = procvrttab[MPI::num_processes()-1] + proccnttab[MPI::num_processes()-1];

  // Perform sanity check
  for (uint i = 1; i <= MPI::process_number(); ++i)
  {
    cout << "Testing procvrttab: " << proccnttab[i-1] << "  " << procvrttab[i-1] << "  " << procvrttab[i] << endl;
    assert( procvrttab[i] >= (procvrttab[i-1] + proccnttab[i-1]) );
  }

}
//-----------------------------------------------------------------------------
