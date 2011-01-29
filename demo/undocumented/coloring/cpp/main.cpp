// Copyright (C) 2010 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2010.
//
// First added:  2010-11-16
// Last changed: 2010-11-17
//
// This demo colors the cells of a mesh such that cells with the same
// color are not neighbors. 'Neighbors' can be in the sense of shared
// vertices, edges or facets.

#include <dolfin.h>

using namespace dolfin;

int main()
{
  // Create mesh
  UnitCube mesh(24, 24, 24);

  // Compute vertex-based coloring
  const MeshFunction<dolfin::uint>& colors_vertex = mesh.color("vertex");
  plot(colors_vertex, "Vertex-based cell coloring");

  // Compute edge-based coloring
  const MeshFunction<dolfin::uint>& colors_edge = mesh.color("edge");
  plot(colors_edge, "Edge-based cell coloring");

  // Compute facet-based coloring
  const MeshFunction<dolfin::uint>& colors_facet = mesh.color("facet");
  plot(colors_facet, "Facet-based cell coloring");

  // Compute facet-based coloring with distance 2
  std::vector<dolfin::uint> coloring_type;
  coloring_type.push_back(mesh.topology().dim());
  coloring_type.push_back(mesh.topology().dim()-1);
  coloring_type.push_back(mesh.topology().dim());
  coloring_type.push_back(mesh.topology().dim()-1);
  coloring_type.push_back(mesh.topology().dim());
  const MeshFunction<dolfin::uint>& colors_vertex_2 = mesh.color(coloring_type);
  plot(colors_vertex_2, "Facet-based cell coloring with distance 2");

  return 0;
}
