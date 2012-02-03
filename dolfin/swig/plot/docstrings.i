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
// Last changed: 2012-02-03

// Autogenerated docstrings file, extracted from the DOLFIN source C++ files.

// Documentation extracted from: (module=plot, header=plot.h)
%feature("docstring")  dolfin::plot "
**Overloaded versions**

* plot\ (v, title=\"Function\", mode=\"auto\")

  Simple built-in plot commands for plotting functions and meshes.
  For plotting to work, PyDOLFIN and Viper must be installed.
  Plot function

* plot\ (v, mesh, title=\"Expression\", mode=\"auto\")

  Plot function

* plot\ (mesh, title=\"Mesh\")

  Plot mesh

* plot\ (f, MeshFunction<uint>\")

  Plot mesh function

* plot\ (f, title=\"MeshFunction<double>\")

  Plot mesh function

* plot\ (f, title=\"MeshFunction<bool>\")

  Plot mesh function
";

// Documentation extracted from: (module=plot, header=FunctionPlotData.h)
%feature("docstring")  dolfin::FunctionPlotData "
This class is used for communicating plot data for functions
to and from (XML) files. It is used by DOLFIN for plotting
Function objects. The data is stored as a mesh and a vector
of interpolated vertex values.
";

%feature("docstring")  dolfin::FunctionPlotData::FunctionPlotData "
**Overloaded versions**

* FunctionPlotData\ (v, mesh)

  Create plot data for given function

* FunctionPlotData\ ()

  Create empty data to be read from file
";

%feature("docstring")  dolfin::FunctionPlotData::vertex_values "
Return vertex values
";

