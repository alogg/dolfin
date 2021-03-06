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

// Autogenerated docstrings file, extracted from the DOLFIN source C++ files.

// Documentation extracted from: (module=quadrature, header=BarycenterQuadrature.h)
%feature("docstring")  dolfin::BarycenterQuadrature "
This class computes the barycenter of an arbitrary polyhedron or
polygon in 3D and therefore allows for barycenter quadrature on
complex polyhedrons. Note: barycenter quadrature is exact for
polynom deg <= 1.
";

%feature("docstring")  dolfin::BarycenterQuadrature::BarycenterQuadrature "
Create barycenter quadrature rule for given CGAL polyhedron
";

%feature("docstring")  dolfin::BarycenterQuadrature::points "
Return points
";

%feature("docstring")  dolfin::BarycenterQuadrature::weights "
Return weights
";

%feature("docstring")  dolfin::BarycenterQuadrature::size "
Return number of quadrature points/weights
";

%feature("docstring")  dolfin::BarycenterQuadrature::compute_quadrature "
Computes barycenter and weight.
";

