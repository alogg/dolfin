// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.4.1.

#ifndef __FFC_09_H
#define __FFC_09_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_09_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_09_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_09_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::tetrahedron;
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
    return 1;
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
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute constants
    const double C0 = element_coordinates[3][0] + element_coordinates[2][0] \
                    + element_coordinates[1][0] - element_coordinates[0][0];
    const double C1 = element_coordinates[3][1] + element_coordinates[2][1] \
                    + element_coordinates[1][1] - element_coordinates[0][1];
    const double C2 = element_coordinates[3][2] + element_coordinates[2][2] \
                    + element_coordinates[1][2] - element_coordinates[0][2];
    
    // Get coordinates and map to the reference (FIAT) element
    double x = coordinates[0];
    double y = coordinates[1];
    double z = coordinates[2];
    
    x = (2.0*d00*x + 2.0*d10*y + 2.0*d20*z - d00*C0 - d10*C1 - d20*C2) / detJ;
    y = (2.0*d01*x + 2.0*d11*y + 2.0*d21*z - d01*C0 - d11*C1 - d21*C2) / detJ;
    z = (2.0*d02*x + 2.0*d12*y + 2.0*d22*z - d02*C0 - d12*C1 - d22*C2) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * (1.0 + x)/(y + z) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * (1.0 + y)/(1.0 - z) - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0;
  }

  /// Evaluate order n derivatives of basis function i at given point in cell
  virtual void evaluate_basis_derivatives(unsigned int i,
                                          unsigned int n,
                                          double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const
  {
    // Extract vertex coordinates
    const double * const * element_coordinates = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = element_coordinates[1][0] - element_coordinates[0][0];
    const double J_01 = element_coordinates[2][0] - element_coordinates[0][0];
    const double J_02 = element_coordinates[3][0] - element_coordinates[0][0];
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    const double J_12 = element_coordinates[3][1] - element_coordinates[0][1];
    const double J_20 = element_coordinates[1][2] - element_coordinates[0][2];
    const double J_21 = element_coordinates[2][2] - element_coordinates[0][2];
    const double J_22 = element_coordinates[3][2] - element_coordinates[0][2];
      
    // Compute sub determinants
    const double d00 = J_11*J_22 - J_12*J_21;
    const double d01 = J_12*J_20 - J_10*J_22;
    const double d02 = J_10*J_21 - J_11*J_20;
    
    const double d10 = J_02*J_21 - J_01*J_22;
    const double d11 = J_00*J_22 - J_02*J_20;
    const double d12 = J_01*J_20 - J_00*J_21;
    
    const double d20 = J_01*J_12 - J_02*J_11;
    const double d21 = J_02*J_10 - J_00*J_12;
    const double d22 = J_00*J_11 - J_01*J_10;
      
    // Compute determinant of Jacobian
    double detJ = J_00*d00 + J_10*d10 + J_20*d20;
    
    // Compute constants
    const double C0 = element_coordinates[3][0] + element_coordinates[2][0] \
                    + element_coordinates[1][0] - element_coordinates[0][0];
    const double C1 = element_coordinates[3][1] + element_coordinates[2][1] \
                    + element_coordinates[1][1] - element_coordinates[0][1];
    const double C2 = element_coordinates[3][2] + element_coordinates[2][2] \
                    + element_coordinates[1][2] - element_coordinates[0][2];
    
    // Get coordinates and map to the reference (FIAT) element
    double x = coordinates[0];
    double y = coordinates[1];
    double z = coordinates[2];
    
    x = (2.0*d00*x + 2.0*d10*y + 2.0*d20*z - d00*C0 - d10*C1 - d20*C2) / detJ;
    y = (2.0*d01*x + 2.0*d11*y + 2.0*d21*z - d01*C0 - d11*C1 - d21*C2) / detJ;
    z = (2.0*d02*x + 2.0*d12*y + 2.0*d22*z - d02*C0 - d12*C1 - d22*C2) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * (1.0 + x)/(y + z) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * (1.0 + y)/(1.0 - z) - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 3;
    
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      combinations[j] = new unsigned int [n];
      for (unsigned int k = 0; k < n; k++)
        combinations[j][k] = 0;
    }
        
    // Generate combinations of derivatives
    for (unsigned int row = 1; row < num_derivatives; row++)
    {
      for (unsigned int num = 0; num < row; num++)
      {
        for (unsigned int col = n-1; col+1 > 0; col--)
        {
          if (combinations[row][col] + 1 > 2)
            combinations[row][col] = 0;
          else
          {
            combinations[row][col] += 1;
            break;
          }
        }
      }
    }
    
    // Compute inverse of Jacobian, components are scaled because of difference in FFC/FIAT reference elements
    const double Jinv[3][3] ={{2*d00 / detJ, 2*d10 / detJ, 2*d20 / detJ}, {2*d01 / detJ, 2*d11 / detJ, 2*d21 / detJ}, {2*d02 / detJ, 2*d12 / detJ, 2*d22 / detJ}};
    
    // Declare transformation matrix
    // Declare pointer to two dimensional array and initialise
    double **transform = new double *[num_derivatives];
        
    for (unsigned int j = 0; j < num_derivatives; j++)
    {
      transform[j] = new double [num_derivatives];
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    }
    
    // Construct transformation matrix
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }
    
    // Reset values
    for (unsigned int j = 0; j < 1*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_z_0 = 1;
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    
    // Table(s) of coefficients
    const static double coefficients0[1][1] = \
    {{1.15470053837925}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[1][1] = \
    {{0}};
    
    const static double dmats1[1][1] = \
    {{0}};
    
    const static double dmats2[1][1] = \
    {{0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0;
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
      }
    }
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives
    delete [] combinations;
    
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    double values[1];
    double coordinates[3];
    
    // Nodal coordinates on reference cell
    static double X[1][3] = {{0.25, 0.25, 0.25}};
    
    // Components for each dof
    static unsigned int components[1] = {0};
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0] - X[i][1] - X[i][2];
    const double w1 = X[i][0];
    const double w2 = X[i][1];
    const double w3 = X[i][2];
    
    // Compute affine mapping x = F(X)
    coordinates[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    coordinates[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    coordinates[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
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
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[0];
    vertex_values[2] = dof_values[0];
    vertex_values[3] = dof_values[0];
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new ffc_09_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_09_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  ffc_09_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~ffc_09_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Discontinuous Lagrange finite element of degree 0 on a tetrahedron";
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
      return false;
      break;
    case 3:
      return true;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = m.num_entities[3];
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
    return 1;
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
    dofs[0] = c.entity_indices[3][0];
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
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
    case 3:
      
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** coordinates,
                                    const ufc::cell& c) const
  {
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new ffc_09_dof_map_0();
  }

};

#endif
