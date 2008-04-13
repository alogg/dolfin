// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.4.4.

#ifndef __FFC_24_H
#define __FFC_24_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_24_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_24_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_24_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Brezzi-Douglas-Marini finite element of degree 1 on a triangle";
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
    return 1;
  }

  /// Return the dimension of the value space for axis i
  virtual unsigned int value_dimension(unsigned int i) const
  {
    return 2;
  }

  /// Evaluate basis function i at given point in cell
  virtual void evaluate_basis(unsigned int i,
                              double* values,
                              const double* coordinates,
                              const ufc::cell& c) const
  {
    throw std::runtime_error("// Function evaluate_basis not generated (compiled with -fno-evaluate_basis)");
  }

  /// Evaluate all basis functions at given point in cell
  virtual void evaluate_basis_all(double* values,
                                  const double* coordinates,
                                  const ufc::cell& c) const
  {
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
  }

  /// Evaluate order n derivatives of basis function i at given point in cell
  virtual void evaluate_basis_derivatives(unsigned int i,
                                          unsigned int n,
                                          double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const
  {
    throw std::runtime_error("// Function evaluate_basis_derivatives not generated (compiled with -fno-evaluate_basis_derivatives)");
  }

  /// Evaluate order n derivatives of all basis functions at given point in cell
  virtual void evaluate_basis_derivatives_all(unsigned int n,
                                              double* values,
                                              const double* coordinates,
                                              const ufc::cell& c) const
  {
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    // The reference points, direction and weights:
    const static double X[6][1][2] = {{{0.666666666666667, 0.333333333333333}}, {{0.333333333333333, 0.666666666666667}}, {{0, 0.333333333333333}}, {{0, 0.666666666666667}}, {{0.333333333333333, 0}}, {{0.666666666666667, 0}}};
    const static double W[6][1] = {{1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[6][1][2] = {{{1, 1}}, {{1, 1}}, {{1, 0}}, {{1, 0}}, {{0, -1}}, {{0, -1}}};
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
      
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
      
    // Compute inverse of Jacobian
    const double Jinv_00 =  J_11 / detJ;
    const double Jinv_01 = -J_01 / detJ;
    const double Jinv_10 = -J_10 / detJ;
    const double Jinv_11 =  J_00 / detJ;
    
    double copyofvalues[2];
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    
    // Compute affine mapping y = F(X)
    double y[2];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
    // Evaluate function at physical points
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Copy old values:
    copyofvalues[0] = values[0];
    copyofvalues[1] = values[1];
    // Do the inverse of div piola 
    values[0] = detJ*(Jinv_00*copyofvalues[0]+Jinv_01*copyofvalues[1]);
    values[1] = detJ*(Jinv_10*copyofvalues[0]+Jinv_11*copyofvalues[1]);
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
  }

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const ufc::cell& c) const
  {
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
      
    // Compute determinant of Jacobian
    double detJ = J_00*J_11 - J_01*J_10;
      
    // Compute inverse of Jacobian
    // Evaluate at vertices and use Piola mapping
    vertex_values[0] = (1.0/detJ)*(dof_values[2]*2*J_00 + dof_values[3]*J_00 + dof_values[4]*(-2*J_01) + dof_values[5]*J_01);
    vertex_values[1] = (1.0/detJ)*(dof_values[0]*2*J_00 + dof_values[1]*J_00 + dof_values[4]*(J_00 + J_01) + dof_values[5]*(2*J_00 - 2*J_01));
    vertex_values[2] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*2*J_01 + dof_values[2]*(J_00 + J_01) + dof_values[3]*(2*J_00 - 2*J_01));
    vertex_values[3] = (1.0/detJ)*(dof_values[2]*2*J_10 + dof_values[3]*J_10 + dof_values[4]*(-2*J_11) + dof_values[5]*J_11);
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*2*J_10 + dof_values[1]*J_10 + dof_values[4]*(J_10 + J_11) + dof_values[5]*(2*J_10 - 2*J_11));
    vertex_values[5] = (1.0/detJ)*(dof_values[0]*J_11 + dof_values[1]*2*J_11 + dof_values[2]*(J_10 + J_11) + dof_values[3]*(2*J_10 - 2*J_11));
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new ffc_24_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_24_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  ffc_24_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~ffc_24_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Brezzi-Douglas-Marini finite element of degree 1 on a triangle";
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
      return true;
      break;
    case 2:
      return false;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = 2*m.num_entities[1];
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

  // Return the geometric dimension of the coordinates this dof map provides
  virtual unsigned int geometric_dimension() const
  {
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 2;
  }

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual unsigned int num_entity_dofs(unsigned int d) const
  {
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = 2*c.entity_indices[1][0];
    dofs[1] = 2*c.entity_indices[1][0] + 1;
    dofs[2] = 2*c.entity_indices[1][1];
    dofs[3] = 2*c.entity_indices[1][1] + 1;
    dofs[4] = 2*c.entity_indices[1][2];
    dofs[5] = 2*c.entity_indices[1][2] + 1;
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   unsigned int facet) const
  {
    switch ( facet )
    {
    case 0:
      dofs[0] = 0;
      dofs[1] = 1;
      break;
    case 1:
      dofs[0] = 2;
      dofs[1] = 3;
      break;
    case 2:
      dofs[0] = 4;
      dofs[1] = 5;
      break;
    }
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(unsigned int* dofs,
                                    unsigned int d, unsigned int i) const
  {
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** coordinates,
                                    const ufc::cell& c) const
  {
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.666666666666667*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[0][1] = 0.666666666666667*x[1][1] + 0.333333333333333*x[2][1];
    coordinates[1][0] = 0.333333333333333*x[1][0] + 0.666666666666667*x[2][0];
    coordinates[1][1] = 0.333333333333333*x[1][1] + 0.666666666666667*x[2][1];
    coordinates[2][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[2][0];
    coordinates[2][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[2][1];
    coordinates[3][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[2][0];
    coordinates[3][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[2][1];
    coordinates[4][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[1][0];
    coordinates[4][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[1][1];
    coordinates[5][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[1][0];
    coordinates[5][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[1][1];
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new ffc_24_dof_map_0();
  }

};

#endif
