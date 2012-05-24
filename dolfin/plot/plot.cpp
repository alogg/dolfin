// Copyright (C) 2007-2009 Anders Logg
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
// Modified by Joachim Berdal Haga, 2008.
// Modified by Garth N. Wells, 2008.
// Modified by Fredrik Valdmanis, 2012.
// Modified by Benjamin Kehlet, 2012
//
// First added:  2007-05-02
// Last changed: 2012-05-22

#include <cstdlib>
#include <sstream>

#include <dolfin/common/MPI.h>
#include <dolfin/common/utils.h>
#include <dolfin/parameter/GlobalParameters.h>
#include <dolfin/io/File.h>
#include <dolfin/log/log.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/Expression.h>
#include "VTKPlotter.h"
#include "FunctionPlotData.h"
#include "plot.h"

using namespace dolfin;

// Template function for plotting objects
template <typename T>
void plot_object(const T& t, std::string title, std::string mode)
{
  info("Plotting %s (%s).",
          t.name().c_str(), t.label().c_str());

  // Don't plot when running in parallel
  if (dolfin::MPI::num_processes() > 1)
  {
    warning("Plotting disabled when running in parallel; see \
https://bugs.launchpad.net/dolfin/+bug/427534");
    return;
  }

  // Get filename prefix
  std::string prefix = parameters["plot_filename_prefix"];

  // Modify prefix and title when running in parallel
  if (dolfin::MPI::num_processes() > 1)
  {
    const dolfin::uint p = dolfin::MPI::process_number();
    prefix += std::string("_p") + to_string(p);
    title += " (process " + to_string(p) + ")";
  }

  // Save to file
  std::string filename = prefix + std::string(".xml");
  File file(filename);
  file << t;

#ifdef HAS_VTK
  VTKPlotter plotter(t);
  plotter.parameters["title"] = title;
  plotter.parameters["mode"] = mode;
  plotter.plot();
#else
  dolfin_error("plot.cpp",
               "plot object",
	       "Plotting disbled. Dolfin has been compiled without VTK support");
#endif
}
//-----------------------------------------------------------------------------
void dolfin::plot(const Function& v,
                  std::string title, std::string mode)
{
  // Duplicate test here since FunctionPlotData may fail in parallel
  // as it does for the eigenvalue demo when vector is not initialized
  // correctly.
  if (dolfin::MPI::num_processes() > 1)
  {
    warning("Plotting disabled when running in parallel; see \
https://bugs.launchpad.net/dolfin/+bug/427534");
    return;
  }

  dolfin_assert(v.function_space()->mesh());
  FunctionPlotData w(v, *v.function_space()->mesh());
  plot_object(w, title, mode);
}
//-----------------------------------------------------------------------------
void dolfin::plot(const Expression& v, const Mesh& mesh,
                  std::string title, std::string mode)
{
  FunctionPlotData w(v, mesh);
  plot_object(w, title, mode);
}
//-----------------------------------------------------------------------------
void dolfin::plot(const Mesh& mesh,
                  std::string title)
{
  plot_object(mesh, title, "auto");
}
//-----------------------------------------------------------------------------
void dolfin::plot(const MeshFunction<uint>& f,
                  std::string title)
{
  plot_object(f, title, "auto");
}
//-----------------------------------------------------------------------------
void dolfin::plot(const MeshFunction<double>& f,
                  std::string title)
{
  plot_object(f, title, "auto");
}
//-----------------------------------------------------------------------------
void dolfin::plot(const MeshFunction<bool>& f,
                  std::string title)
{
  plot_object(f, title, "auto");
}
//-----------------------------------------------------------------------------
