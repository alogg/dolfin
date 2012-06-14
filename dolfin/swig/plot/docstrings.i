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

// Documentation extracted from: (module=plot, header=plot.h)
%feature("docstring")  dolfin::plot "
**Overloaded versions**

* plot\ (v, title=\"Function\", mode=\"auto\")

  Plot function

* plot\ (v, mesh, title=\"Expression\", mode=\"auto\")

  Plot function

* plot\ (mesh, title=\"Mesh\")

  Plot mesh

* plot\ (bc, B.C.\")

  Plot Dirichlet B.C.

* plot\ (f, MeshFunction<uint>\")

  Plot mesh function

* plot\ (f, title=\"MeshFunction<double>\")

  Plot mesh function

* plot\ (f, title=\"MeshFunction<bool>\")

  Plot mesh function
";

%feature("docstring")  dolfin::interactive "
Make the current plot interactive
";

// Documentation extracted from: (module=plot, header=PlotableExpression.h)
%feature("docstring")  dolfin::PlotableExpression "
A light wrapper class to hold an expression to plot, along with the mesh
to plot it on. Allows for clean, templated plotter code in plot.cpp
";

%feature("docstring")  dolfin::PlotableExpression::PlotableExpression "
Create plotable expression object
";

%feature("docstring")  dolfin::PlotableExpression::id "
Return unique ID of the expression
";

%feature("docstring")  dolfin::PlotableExpression::expression "
Get the expression
";

%feature("docstring")  dolfin::PlotableExpression::mesh "
Get the mesh
";

// Documentation extracted from: (module=plot, header=VTKPlotter.h)
%feature("docstring")  dolfin::VTKPlotter "
This class enables visualization of various DOLFIN entities.
It supports visualization of meshes, functions, expressions, boundary
conditions and mesh functions.
The plotter has several parameters that the user can set and adjust to
affect the appearance and behavior of the plot.

A plotter can be created and used in the following way:

  Mesh mesh = ...;
  VTKPlotter plotter(mesh);
  plotter.plot();

Parameters can be adjusted at any time and will take effect on the next
call to the plot() method. The following parameters exist:

============= ============ =============== =================================
 Name          Value type   Default value              Description
============= ============ =============== =================================
 mode           String        \"auto\"        For vector valued functions,
                                            this parameter may be set to
                                            \"warp\" to enable vector warping
                                            visualization
 interactive    Boolean         True        Enable/disable interactive mode
                                            for the rendering window.
                                            For repeated plots of the same
                                            object (animated plots), this
                                            parameter must be set to false
 wireframe      Boolean     True for        Enable/disable wireframe
                            meshes, else    rendering of the object
                            false
 title          String      Inherited       The title of the rendering
                            from the        window
                            name/label of
                            the object
 scale          Double      1.0             Adjusts the scaling of the
                                            warping and glyphs
 scalarbar      Boolean     False for       Hide/show the colormapping bar
                            meshes, else
                            true
============= ============ =============== =================================

The default visualization mode for the different plot types are as follows:

=========================  ============================ ===================
 Plot type                  Default visualization mode   Alternatives
=========================  ============================ ===================
 Meshes                     Wireframe rendering           None
 2D scalar functions        Scalar warping                None
 3D scalar functions        Color mapping                 None
 2D/3D vector functions     Glyphs (vector arrows)        Vector warping
=========================  ============================ ===================

Expressions and boundary conditions are also visualized according to the
above table.
";

%feature("docstring")  dolfin::VTKPlotter::VTKPlotter "
**Overloaded versions**

* VTKPlotter\ (plotter)

  Copy constructor

* VTKPlotter\ (mesh)

  Create plotter for a mesh

* VTKPlotter\ (function)

  Create plotter for a function

* VTKPlotter\ (plotable)

  Create plotter for an expression

* VTKPlotter\ (expression, mesh)

  Create plotter for an expression

* VTKPlotter\ (bc)

  Create plotter for Dirichlet B.C.

* VTKPlotter\ (mesh_function)

  Create plotter for an integer valued mesh function

* VTKPlotter\ (mesh_function)

  Create plotter for a double valued mesh function

* VTKPlotter\ (mesh_function)

  Create plotter for a boolean valued mesh function
";

%feature("docstring")  dolfin::VTKPlotter::operator= "
Assignment operator
";

%feature("docstring")  dolfin::VTKPlotter::plot "
Plot the object
";

%feature("docstring")  dolfin::VTKPlotter::default_parameters "
Default parameter values
";

%feature("docstring")  dolfin::VTKPlotter::default_mesh_parameters "
Default parameter values for mesh plotting
";

%feature("docstring")  dolfin::VTKPlotter::interactive "
Make the current plot interactive
";

%feature("docstring")  dolfin::VTKPlotter::id "
Return unique ID of the object to plot
";

