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

// Documentation extracted from: (module=generation, header=Interval.h)
%feature("docstring")  dolfin::Interval "
Interval mesh of the 1D line [a,b].  Given the number of cells
(nx) in the axial direction, the total number of intervals will
be nx and the total number of vertices will be (nx + 1).
";

%feature("docstring")  dolfin::Interval::Interval "
Constructor

*Arguments*
    nx (int)
        The number of cells.
    a (float)
        The minimum point (inclusive).
    b (float)
        The maximum point (inclusive).

*Example*
    .. note::
    
        No example code available for this function.
";

// Documentation extracted from: (module=generation, header=PolygonalMeshGenerator.h)
%feature("docstring")  dolfin::PolygonalMeshGenerator "
Polygonal mesh generator that uses CGAL
";

%feature("docstring")  dolfin::PolygonalMeshGenerator::generate "
**Overloaded versions**

* generate\ (mesh, vertices, cell_size)

  Generate mesh of a polygonal domain described by domain vertices

* generate\ (mesh, polygon, cell_size)

  Generate mesh of a domain described by a CGAL polygon
";

// Documentation extracted from: (module=generation, header=PolyhedralMeshGenerator.h)
%feature("docstring")  dolfin::PolyhedralMeshGenerator "
Polyhedral mesh generator that uses CGAL
";

%feature("docstring")  dolfin::PolyhedralMeshGenerator::generate "
**Overloaded versions**

* generate\ (mesh, off_file, cell_size, detect_sharp_features=true)

  Create mesh from Object File Format (.off) file

* generate\ (mesh, vertices, facets, cell_size, detect_sharp_features=true)

  Create mesh from a collection of facets
";

%feature("docstring")  dolfin::PolyhedralMeshGenerator::generate_surface_mesh "
Create a surface mesh from Object File Format (.off) file
";

%feature("docstring")  dolfin::PolyhedralMeshGenerator::cgal_generate "
Create mesh from a CGAL polyhedron
";

%feature("docstring")  dolfin::PolyhedralMeshGenerator::cgal_generate_surface_mesh "
Create surface mesh from a CGAL polyhedron
";

// Documentation extracted from: (module=generation, header=Triangulate.h)
%feature("docstring")  dolfin::Triangulate "
Create mesh from a triangulation of points
";

%feature("docstring")  dolfin::Triangulate::triangulate "
Create mesh from a triangulation of points
";

// Documentation extracted from: (module=generation, header=UnitTetrahedron.h)
%feature("docstring")  dolfin::UnitTetrahedron "
A mesh consisting of a single tetrahedron with vertices at

  (0, 0, 0)
  (1, 0, 0)
  (0, 1, 0)
  (0, 0, 1)

This class is useful for testing.
";

%feature("docstring")  dolfin::UnitTetrahedron::UnitTetrahedron "
Create mesh of unit tetrahedron
";

// Documentation extracted from: (module=generation, header=UnitCube.h)
%feature("docstring")  dolfin::UnitCube "
Tetrahedral mesh of the 3D unit cube [0,1] x [0,1] x [0,1].
Given the number of cells (nx, ny, nz) in each direction,
the total number of tetrahedra will be 6*nx*ny*nz and the
total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).
";

%feature("docstring")  dolfin::UnitCube::UnitCube "
Create a uniform finite element :py:class:`Mesh` over the unit cube
[0,1] x [0,1] x [0,1].

*Arguments*
    nx (int)
        Number of cells in :math:`x` direction.
    ny (int)
        Number of cells in :math:`y` direction.
    nz (int)
        Number of cells in :math:`z` direction.

*Example*
    .. note::
    
        No example code available for this function.
";

// Documentation extracted from: (module=generation, header=UnitInterval.h)
%feature("docstring")  dolfin::UnitInterval "
A mesh of the unit interval (0, 1) with a given number of cells
(nx) in the axial direction. The total number of intervals will
be nx and the total number of vertices will be (nx + 1).
";

%feature("docstring")  dolfin::UnitInterval::UnitInterval "
Create mesh of unit interval
";

// Documentation extracted from: (module=generation, header=UnitTriangle.h)
%feature("docstring")  dolfin::UnitTriangle "
A mesh consisting of a single triangle with vertices at

  (0, 0)
  (1, 0)
  (0, 1)

This class is useful for testing.
";

%feature("docstring")  dolfin::UnitTriangle::UnitTriangle "
Create mesh of unit triangle
";

// Documentation extracted from: (module=generation, header=UnitSquare.h)
%feature("docstring")  dolfin::UnitSquare "
Triangular mesh of the 2D unit square [0,1] x [0,1].
Given the number of cells (nx, ny) in each direction,
the total number of triangles will be 2*nx*ny and the
total number of vertices will be (nx + 1)*(ny + 1).

std::string diagonal (\"left\", \"right\", \"right/left\", \"left/right\",
or \"crossed\") indicates the direction of the diagonals.
";

%feature("docstring")  dolfin::UnitSquare::UnitSquare "
Create a uniform finite element :py:class:`Mesh` over the unit square
[0,1] x [0,1].

*Arguments*
    nx (int)
        Number of cells in horizontal direction.
    ny (int)
        Number of cells in vertical direction.
    diagonal (str)
        Optional argument: A std::string indicating
        the direction of the diagonals.

*Example*
    .. note::
    
        No example code available for this function.
";

// Documentation extracted from: (module=generation, header=UnitCircle.h)
%feature("docstring")  dolfin::UnitCircle "
Tetrahedral mesh of the unit circle.
";

%feature("docstring")  dolfin::UnitCircle::UnitCircle "
Create a uniform finite element :py:class:`Mesh` over the unit circle.

*Arguments*
    n (int)
        Resolution of the mesh.
    diagonal (str)
        Optional argument: A std::string indicating
        the direction of the diagonals.
    transformation (str)
        Optional argument: A std::string indicating
        the type of transformation used.
";

// Documentation extracted from: (module=generation, header=Box.h)
%feature("docstring")  dolfin::Box "
Tetrahedral mesh of the 3D rectangular prism [x0, x1] x [y0, y1]
x [z0, z1].  Given the number of cells (nx, ny, nz) in each
direction, the total number of tetrahedra will be 6*nx*ny*nz and
the total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).
";

%feature("docstring")  dolfin::Box::Box "
Create a uniform finite element :py:class:`Mesh` over the rectangular prism
[x0, x1] x [y0, y1] x [z0, z1].

*Arguments*
    x0 (float)
        :math:`x`-min.
    y0 (float)
        :math:`y`-min.
    z0 (float)
        :math:`z`-min.
    x1 (float)
        :math:`x`-max.
    y1 (float)
        :math:`y`-max.
    z1 (float)
        :math:`z`-max.
    xn (float)
        Number of cells in :math:`x`-direction.
    yn (float)
        Number of cells in :math:`y`-direction.
    zn (float)
        Number of cells in :math:`z`-direction.

*Example*
    .. note::
    
        No example code available for this function.
";

// Documentation extracted from: (module=generation, header=Rectangle.h)
%feature("docstring")  dolfin::Rectangle "
Triangular mesh of the 2D rectangle (x0, y0) x (x1, y1).
Given the number of cells (nx, ny) in each direction,
the total number of triangles will be 2*nx*ny and the
total number of vertices will be (nx + 1)*(ny + 1).

*Arguments*
    x0 (float)
        :math:`x`-min.
    y0 (float)
        :math:`y`-min.
    x1 (float)
        :math:`x`-max.
    y1 (float)
        :math:`y`-max.
    xn (float)
        Number of cells in :math:`x`-direction.
    yn (float)
        Number of cells in :math:`y`-direction.
    diagonal (str)
        Direction of diagonals: \"left\", \"right\", \"left/right\", \"crossed\"

*Example*
    .. note::
    
        No example code available for this function.
";

// Documentation extracted from: (module=generation, header=UnitSphere.h)
%feature("docstring")  dolfin::UnitSphere "
Tetrahedral mesh of the unit sphere.
";

%feature("docstring")  dolfin::UnitSphere::UnitSphere "
WARNING:

The UnitSphere class is broken and should not be used for computations.
It generates meshes of very bad quality (very thin tetrahedra).

Create a uniform finite element :py:class:`Mesh` over the unit sphere.

*Arguments*
    n (int)
        Resolution of the mesh.
";

