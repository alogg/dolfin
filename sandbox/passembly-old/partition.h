// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-05-27
// Last changed: 
//
// This file is used for testing parallel assembly

#ifndef __PARTITION_H
#define __PARTITION_H

#include <dolfin.h>

#include <parmetis.h>
extern "C"
{
  #include <metis.h>
}

using namespace dolfin;
//-----------------------------------------------------------------------------
void testMeshPartition(Mesh& mesh, MeshFunction<dolfin::uint>& cell_partition_function,
          MeshFunction<dolfin::uint>& vertex_partition_function, int num_partitions)
{
  int num_cells     = mesh.numCells() ;
  int num_vertices  = mesh.numVertices();
  
  int index_base = 0;  // zero-based indexing
  int edges_cut  = 0;

  int cell_type = 0;
  idxtype* cell_partition   = new int[num_cells];
  idxtype* vertex_partition = new int[num_vertices];
  idxtype* mesh_data = 0;

  // Set cell type and allocate memory for METIS mesh structure
  if(mesh.type().cellType() == CellType::triangle)
  {
    cell_type = 1;
    mesh_data = new int[3*num_cells];
  }
  else if(mesh.type().cellType() == CellType::tetrahedron) 
  {
    cell_type = 2;
    mesh_data = new int[4*num_cells];
  }
  else
    error("Do not know how to partition mesh of this type");
  
  cell_partition_function.init(mesh, mesh.topology().dim());
  vertex_partition_function.init(mesh, 0);

  if(num_partitions > 1)
  {
    // Create mesh structure for METIS
    dolfin::uint i = 0;
    for (CellIterator cell(mesh); !cell.end(); ++cell)
      for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
        mesh_data[i++] = vertex->index();

      // Use METIS to partition mesh
    METIS_PartMeshNodal(&num_cells, &num_vertices, mesh_data, &cell_type, &index_base, 
                        &num_partitions, &edges_cut, cell_partition, vertex_partition);
  
    // Set partition numbers on cells
    i = 0;
    for (CellIterator cell(mesh); !cell.end(); ++cell)
      cell_partition_function.set(cell->index(), cell_partition[i++]);

    // Set partition numbers on vertexes
    i = 0;
    for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
      vertex_partition_function.set(vertex->index(), vertex_partition[i++]);
  }
  else
  {
    // Set partition numbers on cells
    dolfin::uint i = 0;
    for (CellIterator cell(mesh); !cell.end(); ++cell)
      cell_partition_function.set(cell->index(), 0);

    // Set partition numbers on vertexes
    i = 0;
    for (VertexIterator vertex(mesh); !vertex.end(); ++vertex)
      vertex_partition_function.set(vertex->index(), 0);
  
  }
  // Clean up
  delete [] cell_partition;
  delete [] vertex_partition;
  delete [] mesh_data;
}
//-----------------------------------------------------------------------------
#endif
