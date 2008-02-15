// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.4.3.

#ifndef __FFC_10_H
#define __FFC_10_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_10_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_10_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_10_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Discontinuous Lagrange finite element of degree 1 on a tetrahedron";
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
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
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
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = d00*(element_coordinates[0][0] - element_coordinates[2][0] - element_coordinates[3][0]) \
                    + d10*(element_coordinates[0][1] - element_coordinates[2][1] - element_coordinates[3][1]) \
                    + d20*(element_coordinates[0][2] - element_coordinates[2][2] - element_coordinates[3][2]);
    
    const double C1 = d01*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[3][0]) \
                    + d11*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[3][1]) \
                    + d21*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[3][2]);
    
    const double C2 = d02*(element_coordinates[0][0] - element_coordinates[1][0] - element_coordinates[2][0]) \
                    + d12*(element_coordinates[0][1] - element_coordinates[1][1] - element_coordinates[2][1]) \
                    + d22*(element_coordinates[0][2] - element_coordinates[1][2] - element_coordinates[2][2]);
    
    // Get coordinates and map to the UFC reference element
    double x = (C0 + d00*coordinates[0] + d10*coordinates[1] + d20*coordinates[2]) / detJ;
    double y = (C1 + d01*coordinates[0] + d11*coordinates[1] + d21*coordinates[2]) / detJ;
    double z = (C2 + d02*coordinates[0] + d12*coordinates[1] + d22*coordinates[2]) / detJ;
    
    // Map coordinates to the reference cube
    if (std::abs(y + z - 1.0) < 1e-14)
      x = 1.0;
    else
      x = -2.0 * x/(y + z - 1.0) - 1.0;
    if (std::abs(z - 1.0) < 1e-14)
      y = -1.0;
    else
      y = 2.0 * y/(1.0 - z) - 1.0;
    z = 2.0 * z - 1.0;
    
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
    
    // Compute inverse of Jacobian
    const double Jinv[3][3] ={{d00 / detJ, d10 / detJ, d20 / detJ}, {d01 / detJ, d11 / detJ, d21 / detJ}, {d02 / detJ, d12 / detJ, d22 / detJ}};
    
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
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_10_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    
    // Table(s) of coefficients
    const static double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0, 0, 0.223606797749979}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[4][4] = \
    {{0, 0, 0, 0},
    {6.32455532033676, 0, 0, 0},
    {0, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats1[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {5.47722557505166, 0, 0, 0},
    {0, 0, 0, 0}};
    
    const static double dmats2[4][4] = \
    {{0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0},
    {1.82574185835055, 0, 0, 0},
    {5.16397779494322, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3];
        }
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0];
          new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1];
          new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2];
          new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3;
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
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      delete [] combinations[row];
      delete [] transform[row];
    }
    
    delete [] combinations;
    delete [] transform;
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
    const static double X[4][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[4][1] = {{1}, {1}, {1}, {1}};
    const static double D[4][1][1] = {{{1}}, {{1}}, {{1}}, {{1}}};
    
    const double * const * x = c.coordinates;
    double result = 0.0;
    // Iterate over the points:
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0][0] - X[i][0][1] - X[i][0][2];
    const double w1 = X[i][0][0];
    const double w2 = X[i][0][1];
    const double w3 = X[i][0][2];
    
    // Compute affine mapping y = F(X)
    double y[3];
    y[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
    y[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
    y[2] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];
    
    // Evaluate function at physical points
    double values[1];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 1; k++)
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
    // Evaluate at vertices and use affine mapping
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
    return new ffc_10_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_10_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  ffc_10_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~ffc_10_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Discontinuous Lagrange finite element of degree 1 on a tetrahedron";
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
    __global_dimension = 4*m.num_entities[3];
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

  // Return the geometric dimension of the coordinates this dof map provides
  virtual unsigned int geometric_dimension() const
  {
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 0;
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
    dofs[0] = 4*c.entity_indices[3][0];
    dofs[1] = 4*c.entity_indices[3][0] + 1;
    dofs[2] = 4*c.entity_indices[3][0] + 2;
    dofs[3] = 4*c.entity_indices[3][0] + 3;
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
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[0][2] = x[0][2];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[1][2] = x[1][2];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[2][2] = x[2][2];
    coordinates[3][0] = x[3][0];
    coordinates[3][1] = x[3][1];
    coordinates[3][2] = x[3][2];
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new ffc_10_dof_map_0();
  }

};

#endif
