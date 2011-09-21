// Copyright (C) 2008-2009 Niclas Jansson, Ola Skavhaug and Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Garth N. Wells, 2010
// Modified by Kent-Andre Mardal, 2011
//
// First added:  2008-12-01
// Last changed: 2010-04-04

#ifndef __MESH_PARTITIONING_H
#define __MESH_PARTITIONING_H

#include <map>
#include <utility>
#include <vector>
#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include "LocalMeshValueCollection.h"
#include "MeshPartitioning.h"

namespace dolfin
{

  template <class T> class LocalMeshValueCollection;
  class Mesh;
  template <class T> class MeshFunction;
  template <class T> class MeshValueCollection;
  class LocalMeshData;

  /// This class partitions and distributes a mesh based on
  /// partitioned local mesh data. Note that the local mesh data will
  /// also be repartitioned and redistributed during the computation
  /// of the mesh partitioning.
  ///
  /// After partitioning, each process has a local mesh and set of
  /// mesh data that couples the meshes together.
  ///
  /// The following mesh data is created:
  ///
  /// 1. "global entity indices 0" (MeshFunction<uint>)
  ///
  /// This maps each local vertex to its global index.
  ///
  /// 2. "overlap" (std::map<uint, std::vector<uint> >)
  ///
  /// This maps each shared vertex to a list of the processes sharing
  /// the vertex.
  ///
  /// 3. "global entity indices %d" (MeshFunction<uint>)
  ///
  /// After partitioning, the function number_entities() may be called
  /// to create global indices for all entities of a given topological
  /// dimension. These are stored as mesh data (MeshFunction<uint>)
  /// named
  ///
  ///    "global entity indices 1"
  ///    "global entity indices 2"
  ///    etc
  ///
  /// 4. "num global entities" (std::vector<uint>)
  ///
  /// The function number_entities also records the number of global
  /// entities for the dimension of the numbered entities in the array
  /// named "num global entities". This array has size D + 1, where D
  /// is the topological dimension of the mesh. This array is
  /// initially created by the mesh and then contains only the number
  /// entities of dimension 0 (vertices) and dimension D (cells).

  class MeshPartitioning
  {
  public:

    template<class T>
    static void build_distributed_value_collection(MeshValueCollection<T>& values,const Mesh& mesh) {}

   /// Build a partitioned mesh based on local meshes
    static void build_distributed_mesh(Mesh& mesh);

    /// Build a partitioned mesh based on local mesh data
    static void build_distributed_mesh(Mesh& mesh, LocalMeshData& data);

    template<class T>
    static void build_distributed_value_collection(MeshValueCollection<T>& values,
               const LocalMeshValueCollection<T>& local_data, const Mesh& mesh);

    /// Create global entity indices for entities of dimension d
    static void number_entities(const Mesh& mesh, uint d);

  private:

    /// Create a partitioned mesh based on local mesh data
    static void partition(Mesh& mesh, LocalMeshData& data);

    /// Create and attach distributed MeshDomains from local_data
    static void build_mesh_domains(Mesh& mesh, const LocalMeshData& local_data);

    /// Create and attach distributed MeshDomains from local_data
    /// [entry, (cell_index, local_index, value)]
    template<class T>
    static void build_mesh_value_collection(const Mesh& mesh,
         const std::vector<std::pair<std::pair<uint, uint>, T> >& local_value_data,
         MeshValueCollection<T>& mesh_values);

    // Compute and return (number of global entities, process offset)
    static std::pair<uint, uint> compute_num_global_entities(uint num_local_entities,
                                                     uint num_processes,
                                                     uint process_number);

    // Build preliminary 'guess' of shared enties
    static void compute_preliminary_entity_ownership(const std::map<std::vector<uint>, uint>& entities,
          const std::map<uint, std::vector<uint> >& shared_vertices,
          std::map<std::vector<uint>, uint>& owned_entity_indices,
          std::map<std::vector<uint>, uint>& shared_entity_indices,
          std::map<std::vector<uint>, std::vector<uint> >& shared_entity_processes,
          std::map<std::vector<uint>, uint>& ignored_entity_indices,
          std::map<std::vector<uint>, std::vector<uint> >& ignored_entity_processes);

    // Communicate with other processes to finalise entity ownership
    static void compute_final_entity_ownership(std::map<std::vector<uint>, uint>& owned_entity_indices,
          std::map<std::vector<uint>, uint>& shared_entity_indices,
          std::map<std::vector<uint>, std::vector<uint> >& shared_entity_processes,
          std::map<std::vector<uint>, uint>& ignored_entity_indices,
          std::map<std::vector<uint>, std::vector<uint> >& ignored_entity_processes);

    // Distribute cells
    static void distribute_cells(LocalMeshData& data,
                                 const std::vector<uint>& cell_partition);

    // Distribute vertices
    static void distribute_vertices(LocalMeshData& data,
                                    std::map<uint, uint>& glob2loc);

    // Build mesh
    static void build_mesh(Mesh& mesh, const LocalMeshData& data,
                           std::map<uint, uint>& glob2loc);

    // Check if all entity vertices are in overlap
    static bool in_overlap(const std::vector<uint>& entity_vertices,
                           const std::map<uint, std::vector<uint> >& overlap);

    // Mark non-shared mesh entities
    static void mark_nonshared(const std::map<std::vector<uint>, uint>& entities,
               const std::map<std::vector<uint>, uint>& shared_entity_indices,
               const std::map<std::vector<uint>, uint>& ignored_entity_indices,
               MeshFunction<bool>& exterior_facets);
  };

  //---------------------------------------------------------------------------
  template<class T>
  inline void MeshPartitioning::build_distributed_value_collection(MeshValueCollection<T>& values,
             const LocalMeshValueCollection<T>& local_data, const Mesh& mesh)
  {
    error("MeshPartitioning::build_distributed_value_collection not implemented.");
  }
  //---------------------------------------------------------------------------

}

#endif
