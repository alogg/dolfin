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

// Documentation extracted from: (module=ale, header=ALE.h)
%feature("docstring")  dolfin::ALE "
This class provides functionality useful for implementation of
ALE (Arbitrary Lagrangian-Eulerian) methods, in particular
moving the boundary vertices of a mesh and then interpolating
the new coordinates for the interior vertices accordingly.
";

%feature("docstring")  dolfin::ALE::move "
**Overloaded versions**

* move\ (mesh, new_boundary)

  Move coordinates of mesh according to new boundary coordinates

* move\ (mesh0, mesh1)

  Move coordinates of mesh0 according to mesh1 with common global vertices

* move\ (mesh, displacement)

  Move coordinates of mesh according to displacement function
";

