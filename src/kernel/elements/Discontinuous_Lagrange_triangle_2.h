// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.3.5.

#ifndef __DISCONTINUOUS_LAGRANGE_TRIANGLE_2_H
#define __DISCONTINUOUS_LAGRANGE_TRIANGLE_2_H

#include <ufc.h>

/// This class defines the interface for a finite element.

class Discontinuous_Lagrange_triangle_2_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  Discontinuous_Lagrange_triangle_2_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Discontinuous_Lagrange_triangle_2_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Discontinuous Lagrange finite element of degree 2 on a triangle";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::triangle;
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
    return 6;
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
    // Not implemented
    return 0.0;
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[2];
    vertex_values[2] = dof_values[5];
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new Discontinuous_Lagrange_triangle_2_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Discontinuous_Lagrange_triangle_2_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Discontinuous_Lagrange_triangle_2_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Discontinuous_Lagrange_triangle_2_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Discontinuous Lagrange finite element of degree 2 on a triangle";
  }

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const
  {
    switch ( d )
    {
    case 0:
      return false;
      break;
    case 1:
      return false;
      break;
    case 2:
      return true;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = 6*m.num_entities[2];
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
    return 6;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 0;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = 6*c.entity_indices[2][0];
    dofs[1] = 6*c.entity_indices[2][0] + 1;
    dofs[2] = 6*c.entity_indices[2][0] + 2;
    dofs[3] = 6*c.entity_indices[2][0] + 3;
    dofs[4] = 6*c.entity_indices[2][0] + 4;
    dofs[5] = 6*c.entity_indices[2][0] + 5;
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
      
      break;
    case 1:
      
      break;
    case 2:
      
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
