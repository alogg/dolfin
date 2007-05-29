// Copyright (C) 2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-05-29
// Last changed: 
//
// This file is used for testing parallel assembly

#ifndef __PDOF_MAP_H
#define __PDOF_MAP_H

#include <set>
#include <dolfin.h>
#include <dolfin/UFC.h>

namespace dolfin
{
  // This is class re-computes the dof map taking into account the mesh partition

  class PdofMap
  {
  public:

    PdofMap(unsigned int size, unsigned int num_processes) 
        : map(size), _process_dimension(num_processes)
    {
      map.clear();
    };

    ~PdofMap()
    {
      // Do nothing
    };

    void create(unsigned int num_processes, unsigned int this_process, Mesh& mesh, UFC& ufc,
                MeshFunction<unsigned int>& cell_partition_function)
    {
      std::set<unsigned int> set;
      std::pair<std::set<unsigned int>::const_iterator, bool> set_return;

      unsigned int dof_counter = 0;
      _process_dimension.clear();

      for (unsigned int proc = 0; proc < num_processes; ++proc)
      {
        for (CellIterator cell(mesh); !cell.end(); ++cell)
        {
          if (cell_partition_function.get(*cell) != static_cast<unsigned int>(proc))
            continue;

          // Update to current cell
          ufc.update(*cell);
  
          // Tabulate dofs for each dimension
          ufc.dof_maps[0]->tabulate_dofs(ufc.dofs[0], ufc.mesh, ufc.cell);
          ufc.dof_maps[1]->tabulate_dofs(ufc.dofs[1], ufc.mesh, ufc.cell);

          for(unsigned int i=0; i < ufc.dof_maps[0]->local_dimension(); ++i)
          {
            set_return = set.insert( (ufc.dofs[0])[i] );
            if( set_return.second )
            {
              map[ (ufc.dofs[0])[i] ] = dof_counter++;
              _process_dimension[proc]++;
            }
          }
        }
      }
    };

    unsigned int process_dimensions(unsigned int dim, unsigned int process)
      { return _process_dimension[process]; }; 

    void update(unsigned int dofs[], unsigned int ufc_dofs[], unsigned int dimension)
    {
      for(unsigned int i = 0; i < dimension; ++i)
        dofs[i] = ufc_dofs[i];
    };

  private:

    std::vector<unsigned int> map;
    std::vector<unsigned int> _process_dimension;

  };
}
//-----------------------------------------------------------------------------
#endif
