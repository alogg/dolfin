// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.3.5.

#ifndef __LAGRANGE_TETRAHEDRON_1_H
#define __LAGRANGE_TETRAHEDRON_1_H

#include <ufc.h>

/// This class defines the interface for a finite element.

class Lagrange_tetrahedron_1_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  Lagrange_tetrahedron_1_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Lagrange_tetrahedron_1_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Lagrange finite element of degree 1 on a tetrahedron";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::tetrahedron;
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
    return 4;
  }

  /// Return the rank of the value space
  virtual unsigned int value_rank() const
  {
    return 0;
  }

  /// Return the dimension of the value space for axis i
  virtual unsigned int value_dimension(unsigned int i) const
  {
    return 1;
  }

  /// Evaluate basis function i at given point in cell
  virtual void evaluate_basis(unsigned int i,
                              double* values,
                              const double* coordinates,
                              const ufc::cell& c) const
  {
    // Not implemented
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    static double values[1];
    static double coordinates[3];
    
    // Nodal coordinates on reference cell
    static double X[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    
    // Components for each dof
    static unsigned int components[4] = {0, 0, 0, 0};
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0] - X[i][1] - X[i][2];
    const double w1 = X[i][0];
    const double w2 = X[i][1];
    const double w3 = X[i][2];
    
    // Compute affine mapping x = F(X)
    coordinates[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    coordinates[0] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    coordinates[0] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at coordinates
    f.evaluate(values, coordinates, c);
    
    // Pick component for evaluation
    return values[components[i]];
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
    vertex_values[3] = dof_values[3];
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new Lagrange_tetrahedron_1_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Lagrange_tetrahedron_1_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Lagrange_tetrahedron_1_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Lagrange_tetrahedron_1_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Lagrange finite element of degree 1 on a tetrahedron";
  }

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const
  {
    switch ( d )
    {
    case 0:
      return true;
      break;
    case 1:
      return false;
      break;
    case 2:
      return false;
      break;
    case 3:
      return false;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = m.num_entities[0];
    return false;
  }

  /// Initialize dof map for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c)
  {
    // Do nothing
  }

  /// Finish initialization of dof map for cells
  virtual void init_cell_finalize()
  {
    // Do nothing
  }

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const
  {
    return __global_dimension;
  }

  /// Return the dimension of the local finite element function space
  virtual unsigned int local_dimension() const
  {
    return 4;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 3;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   const ufc::mesh& m,
                                   const ufc::cell& c,
                                   unsigned int facet) const
  {
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::mesh& m,
                                    const ufc::cell& c) const
  {
    // Not implemented
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    // Not implemented
    return 0;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    // Not implemented
    return 0;
  }

};

#endif
