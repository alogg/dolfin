// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007
// Modified by Johan Hake 2009
//
// First added:  2007-07-08
// Last changed: 2009-09-16

#include <boost/scoped_array.hpp>
#include <vector>
#include <map>

#include <dolfin/log/log.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/common/constants.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/SubDomain.h>
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include "DofMap.h"
#include "UFCMesh.h"
#include "BoundaryCondition.h"
#include "PeriodicBC.h"

using namespace dolfin;

// Comparison operator for hashing coordinates. Note that two
// coordinates are considered equal if equal to within round-off.
struct lt_coordinate
{
  bool operator() (const std::vector<double>& x, const std::vector<double>& y) const
  {
    unsigned int n = std::max(x.size(), y.size());
    for (unsigned int i = 0; i < n; ++i)
    {
      double xx = 0.0;
      double yy = 0.0;

      if (i < x.size())
        xx = x[i];
      if (i < y.size())
        yy = y[i];

      if (xx < (yy - DOLFIN_EPS))
        return true;
      else if (xx > (yy + DOLFIN_EPS))
        return false;
    }

    return false;
  }
};

// Mapping from coordinates to dof pairs
typedef std::map<std::vector<double>, std::pair<int, int>, lt_coordinate> coordinate_map;
typedef coordinate_map::iterator coordinate_iterator;

//-----------------------------------------------------------------------------
PeriodicBC::PeriodicBC(const FunctionSpace& V,
                       const SubDomain& sub_domain)
  : BoundaryCondition(V), sub_domain(reference_to_no_delete_pointer(sub_domain)),
    num_dof_pairs(0), master_dofs(0), slave_dofs(0), zeros(0)
{
  not_working_in_parallel("Periodic boundary conditions");

  // Build mapping
  rebuild();
}
//-----------------------------------------------------------------------------
PeriodicBC::PeriodicBC(boost::shared_ptr<const FunctionSpace> V,
                       boost::shared_ptr<const SubDomain> sub_domain)
  : BoundaryCondition(V), sub_domain(sub_domain),
    num_dof_pairs(0), master_dofs(0), slave_dofs(0), zeros(0)
{
  not_working_in_parallel("Periodic boundary conditions");

  // Build mapping
  rebuild();
}
//-----------------------------------------------------------------------------
PeriodicBC::~PeriodicBC()
{
  delete [] master_dofs;
  delete [] slave_dofs;
  delete [] zeros;
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericMatrix& A) const
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericVector& b) const
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericMatrix& A, GenericVector& b) const
{
  assert(num_dof_pairs > 0);
  assert(master_dofs);
  assert(slave_dofs);
  assert(zeros);

  cout << "Applying periodic boundary conditions to linear system." << endl;

  // Add slave rows to master rows
  std::vector<uint> columns;
  std::vector<double> values;
  for (uint i = 0; i < num_dof_pairs; ++i)
  {
    // Add slave row to master row in A
    A.getrow(slave_dofs[i], columns, values);
    A.add(&values[0], 1, &master_dofs[i], columns.size(), &columns[0]);

    // Add slave row to master row in b
    b.get(&values[0], 1, &slave_dofs[i]);
    b.add(&values[0], 1, &master_dofs[i]);

    // Apply, needed between calls to get and add
    A.apply();
    b.apply();
  }

  // Zero slave rows and insert 1 on the diagonal
  A.ident(num_dof_pairs, slave_dofs);
  b.set(zeros, num_dof_pairs, slave_dofs);

  // Insert -1 for master dofs in slave rows
  for (uint i = 0; i < num_dof_pairs; ++i)
  {
    const double value = -1;
    A.set(&value, 1, &slave_dofs[i], 1, &master_dofs[i]);
  }

  A.apply();
  b.apply();
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericVector& b, const GenericVector& x) const
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
void PeriodicBC::apply(GenericMatrix& A,
                       GenericVector& b,
                       const GenericVector& x) const
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
void PeriodicBC::rebuild()
{
  // FIXME: Make this work for non-scalar subsystems, like vector-valued
  // FIXME: Lagrange where more than one per element is associated with
  // FIXME: each coordinate. Note that globally there may very well be
  // FIXME: more than one dof per coordinate (for conforming elements).

  cout << "Building mapping between periodic degrees of freedom." << endl;

  // Get mesh and dofmap
  assert(V);
  const Mesh& mesh = V->mesh();
  const DofMap& dofmap = V->dofmap();

  // Get dimensions
  const uint tdim = mesh.topology().dim();
  const uint gdim = mesh.geometry().dim();

  // Set geometric dimension (needed for SWIG interface)
  sub_domain->_geometric_dimension = gdim;

  // Make sure we have the facet - cell connectivity
  mesh.init(tdim - 1, tdim);

  // Create local data for application of boundary conditions
  BoundaryCondition::LocalData data(*V);

  // Arrays used for mapping coordinates
  std::vector<double> xx(gdim);
  boost::scoped_array<double> y(new double[gdim]);

  // Mapping from coordinates to dof pairs
  coordinate_map coordinate_dof_pairs;

  // Iterate over all facets of the mesh (not only the boundary)
  Progress p("Applying periodic boundary conditions", mesh.size(tdim - 1));
  for (FacetIterator facet(mesh); !facet.end(); ++facet)
  {
    // Get cell (there may be two, but pick first) and local facet index
    Cell cell(mesh, facet->entities(tdim)[0]);
    const uint local_facet = cell.index(*facet);

    // Tabulate dofs and coordinates on cell
    dofmap.tabulate_dofs(data.cell_dofs, cell);
    dofmap.tabulate_coordinates(data.coordinates, cell);

    // Tabulate which dofs are on the facet
    dofmap.tabulate_facet_dofs(data.facet_dofs, local_facet);

    // Iterate over facet dofs
    for (uint i = 0; i < dofmap.num_facet_dofs(); ++i)
    {
      // Get dof and coordinate of dof
      const uint local_dof = data.facet_dofs[i];
      const int global_dof = data.cell_dofs[local_dof];
      const double* x = data.coordinates[local_dof];

      // Set y = x
      for (uint j = 0; j < gdim; ++j)
        y[j] = x[j];

      // Map coordinate from H to G
      sub_domain->map(x, y.get());

      // Check if coordinate is inside the domain G or in H
      const bool on_boundary = facet->num_entities(tdim) == 1;
      if (sub_domain->inside(x, on_boundary))
      {
        // Copy coordinate to std::vector
        for (uint j = 0; j < gdim; ++j)
          xx[j] = x[j];

        // Check if coordinate exists from before
        coordinate_iterator it = coordinate_dof_pairs.find(xx);
        if (it != coordinate_dof_pairs.end())
        {
          // Exists from before, so set dof associated with x
          it->second.first = global_dof;
        }
        else
        {
          // Doesn't exist, so create new pair (with illegal second value)
          std::pair<int, int> dofs(global_dof, -1);
          coordinate_dof_pairs[xx] = dofs;
        }
      }
      else if (sub_domain->inside(y.get(), on_boundary))
      {
        // Copy coordinate to std::vector
        for (uint j = 0; j < gdim; ++j)
          xx[j] = y[j];

        // Check if coordinate exists from before
        coordinate_iterator it = coordinate_dof_pairs.find(xx);
        if (it != coordinate_dof_pairs.end())
        {
          // Exists from before, so set dof associated with y
          it->second.second = global_dof;
        }
        else
        {
          // Doesn't exist, so create new pair (with illegal first value)
          std::pair<int, int> dofs(-1, global_dof);
          coordinate_dof_pairs[xx] = dofs;
        }
      }
    }

    p++;
  }

  // Delete old arrays if necessary
  if (master_dofs)
    delete [] master_dofs;
  if (slave_dofs)
    delete [] slave_dofs;
  if (zeros)
    delete [] zeros;

  // Initialize arrays
  num_dof_pairs = coordinate_dof_pairs.size();
  master_dofs = new uint[num_dof_pairs];
  slave_dofs = new uint[num_dof_pairs];
  zeros = new double[num_dof_pairs];
  uint pos = 0;

  // Store master and slave dofs
  for (coordinate_iterator it = coordinate_dof_pairs.begin(); it != coordinate_dof_pairs.end(); ++it)
  {
    // Check dofs
    if (it->second.first == -1 || it->second.second == -1)
    {
      cout << "At coordinate: x =";
      for (uint j = 0; j < gdim; ++j)
        cout << " " << it->first[j];
      cout << endl;
      error("Unable to find a pair of matching dofs for periodic boundary condition.");
    }

    // Store dofs
    master_dofs[pos] = it->second.first;
    slave_dofs[pos] = it->second.second;
    zeros[pos] = 0.0;

    pos++;
  }
}
//-----------------------------------------------------------------------------
