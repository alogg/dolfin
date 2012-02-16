// Copyright (C) 2012 Garth N. Wells
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
// First added:  2012-02-03
// Last changed: 2012-02-16

#ifndef __DOLFIN_POLYHEDRALMESHGENERATOR_H
#define __DOLFIN_POLYHEDRALMESHGENERATOR_H

#ifdef HAS_CGAL

#include <string>
#include <vector>

namespace dolfin
{

  // Forward declarations
  class Mesh;

  /// Polyhedral mesh generator that uses CGAL

  class PolyhedralMeshGenerator
  {
  public:

    /// Create mesh from Object File Format (.off) file
    static void generate(Mesh& mesh, const std::string off_file,
                         double cell_size, bool detect_sharp_features=true);

    /// Create mesh from a collection of facets
    static void generate(Mesh& mesh, const std::vector<Point>& vertices,
                         const std::vector<std::vector<unsigned int> >& facets,
                         double cell_size, bool detect_sharp_features=true);

  private:

    /// Create mesh from a CGAL mesh domain
    template<typename T>
    static void cgal_generate(Mesh& mesh, T& p, double cell_size,
                              bool detect_sharp_features);
  };

}

#endif
#endif
