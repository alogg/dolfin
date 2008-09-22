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
      dolfin_debug("create()");
      std::set<unsigned int> set;
      std::pair<std::set<unsigned int>::const_iterator, bool> set_return;

      unsigned int dof_counter = 0;
      _process_dimension.clear();

      for (unsigned int proc = 0; proc < num_processes; ++proc)
      {
        dolfin_debug1("proc = %d", proc);
        for (CellIterator cell(mesh); !cell.end(); ++cell)
        {
          if (cell_partition_function.get(*cell) != static_cast<unsigned int>(proc))
            continue;

          dolfin_debug1("Updating cell = %d", cell->index());
          // Update to current cell
          ufc.update(*cell);
  
          // Tabulate dofs for each dimension
          dolfin_debug1("Tabulate dofs for cell = %d dim 0", cell->index());
          ufc_dof_map[0].tabulate_dofs(ufc.dofs[0], ufc.mesh, ufc.cell);
          dolfin_debug1("Tabulate dofs for cell = %d dim 1", cell->index());
          ufc_dof_map[1].tabulate_dofs(ufc.dofs[1], ufc.mesh, ufc.cell);

          for(unsigned int i=0; i < ufc_dof_map->local_dimension(); ++i)
          {
            dolfin_debug1("dof = %d", i);
            set_return = set.insert( (ufc.dofs[0])[i] );
            if( set_return.second )
            {
              map[ (ufc.dofs[0])[i] ] = dof_counter++;
              _process_dimension[proc]++;
            }
          }
        }
      }
      for(unsigned int i=0; i<map.size(); ++i)
      {
        dolfin_debug2("map[%d] = %d", i, map[i]);
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

    ufc::dof_map* ufc_dof_map;

  };
}
//-----------------------------------------------------------------------------
#endif
