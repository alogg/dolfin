// Copyright (C) 2011 Anders Logg
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
// First added:  2011-08-29
// Last changed: 2011-08-30

#ifndef __MESH_DOMAINS_H
#define __MESH_DOMAINS_H

#include <vector>
#include <boost/shared_ptr.hpp>

#include <dolfin/common/types.h>
#include "MeshMarkers.h"

namespace dolfin
{

  // Forward declarations
  class Mesh;
  template <class T> class MeshFunction;

  /// The class _MeshDomains_ stores the division of a _Mesh_ into
  /// subdomains. For each topological dimension 0 <= d <= D, where D
  /// is the topological dimension of the _Mesh_, a set of integer
  /// markers are stored for a subset of the entities of dimension d,
  /// indicating for each entity in the subset the number of the
  /// subdomain. It should be noted that the subset does not need to
  /// contain all entities of any given dimension; entities not
  /// contained in the subset are "unmarked".

  class MeshDomains
  {
  public:

    /// Create empty mesh domains
    MeshDomains();

    /// Destructor
    ~MeshDomains();

    /// Return number of marked entities of given dimension
    uint num_marked(uint dim) const;

    /// Clear all data
    void clear();

  private:

    // Initialize mesh function of dimension d from markers
    void init_subdomains(uint d);

    // The mesh
    boost::shared_ptr<Mesh> _mesh;

    // Subdomain markers (input/storage)
    std::vector<std::vector<std::vector<uint> > > _markers;

    std::vector<MeshMarkers<uint> > __markers;

    // Subdomains corresponding to markers
    std::vector<boost::shared_ptr<MeshFunction<uint> > > _subdomains;

  };

}

#endif
