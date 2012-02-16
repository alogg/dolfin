// Auto generated SWIG file for Python interface of DOLFIN
//
// Copyright (C) 2012 Kristian B. Oelgaard
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
// First added:  2012-01-18
// Last changed: 2012-02-14

// Autogenerated docstrings file, extracted from the DOLFIN source C++ files.

// Documentation extracted from: (module=generation, header=PolygonalMeshGenerator.h)
%feature("docstring")  dolfin::PolygonalMeshGenerator "
Polygonal mesh generator that uses CGAL
";

%feature("docstring")  dolfin::PolygonalMeshGenerator::generate "
Generate mesh of a polygonal domain described by domain vertices
";

// Documentation extracted from: (module=generation, header=PolyhedralMeshGenerator.h)
%feature("docstring")  dolfin::PolyhedralMeshGenerator "
Polyhedral mesh generator that uses CGAL
";

%feature("docstring")  dolfin::PolyhedralMeshGenerator::generate "
**Overloaded versions**

* generate\ (mesh, off_file, cell_size)

  Create mesh from Object File Format (.off) file

* generate\ (mesh, vertices, facets, cell_size)

  Create mesh from a collection of facets
";

%feature("docstring")  dolfin::PolyhedralMeshGenerator::cgal_generate "
Create mesh from a CGAL mesh domain
";

// Documentation extracted from: (module=generation, header=Triangulate.h)
%feature("docstring")  dolfin::Triangulate "
Create mesh from a triangulation of points
";

%feature("docstring")  dolfin::Triangulate::triangulate "
Create mesh from a triangulation of points
";

