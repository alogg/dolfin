// This code conforms with the UFC specification version 1.4
// and was automatically generated by FFC version 0.9.2.
//
// This code was generated with the option '-l dolfin' and
// contains DOLFIN-specific wrappers that depend on DOLFIN.
// 
// This code was generated with the following parameters:
// 
//   cache_dir:                      ''
//   convert_exceptions_to_warnings: False
//   cpp_optimize:                   False
//   epsilon:                        1e-14
//   form_postfix:                   True
//   format:                         'dolfin'
//   log_level:                      10
//   log_prefix:                     ''
//   optimize:                       False
//   output_dir:                     '.'
//   precision:                      15
//   quadrature_degree:              'auto'
//   quadrature_rule:                'auto'
//   representation:                 'auto'
//   split:                          False

#ifndef __MASSMATRIX3D_H
#define __MASSMATRIX3D_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>

/// This class defines the interface for a finite element.

class massmatrix3d_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  massmatrix3d_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~massmatrix3d_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "FiniteElement('Lagrange', Cell('tetrahedron', 1, Space(3)), 1)";
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
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_02 = x[3][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    const double J_12 = x[3][1] - x[0][1];
    const double J_20 = x[1][2] - x[0][2];
    const double J_21 = x[2][2] - x[0][2];
    const double J_22 = x[3][2] - x[0][2];
    
    // Compute sub determinants
    const double d_00 = J_11*J_22 - J_12*J_21;
    const double d_01 = J_12*J_20 - J_10*J_22;
    const double d_02 = J_10*J_21 - J_11*J_20;
    const double d_10 = J_02*J_21 - J_01*J_22;
    const double d_11 = J_00*J_22 - J_02*J_20;
    const double d_12 = J_01*J_20 - J_00*J_21;
    const double d_20 = J_01*J_12 - J_02*J_11;
    const double d_21 = J_02*J_10 - J_00*J_12;
    const double d_22 = J_00*J_11 - J_01*J_10;
    
    // Compute determinant of Jacobian
    double detJ = J_00*d_00 + J_10*d_10 + J_20*d_20;
    
    // Compute inverse of Jacobian
    
    // Compute constants
    const double C0 = x[3][0] + x[2][0] + x[1][0] - x[0][0];
    const double C1 = x[3][1] + x[2][1] + x[1][1] - x[0][1];
    const double C2 = x[3][2] + x[2][2] + x[1][2] - x[0][2];
    
    // Get coordinates and map to the reference (FIAT) element
    double X = (d_00*(2.0*coordinates[0] - C0) + d_10*(2.0*coordinates[1] - C1) + d_20*(2.0*coordinates[2] - C2)) / detJ;
    double Y = (d_01*(2.0*coordinates[0] - C0) + d_11*(2.0*coordinates[1] - C1) + d_21*(2.0*coordinates[2] - C2)) / detJ;
    double Z = (d_02*(2.0*coordinates[0] - C0) + d_12*(2.0*coordinates[1] - C1) + d_22*(2.0*coordinates[2] - C2)) / detJ;
    
    
    // Reset values.
    *values = 0.000000000000000;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Array of basisvalues.
    double basisvalues[4] = {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000};
    
    // Declare helper variables.
    unsigned int rr = 0;
    unsigned int ss = 0;
    double tmp0 = 0.500000000000000*(2.000000000000000 + Y + Z + 2.000000000000000*X);
    
    // Compute basisvalues.
    basisvalues[0] = 1.000000000000000;
    basisvalues[1] = tmp0;
    for (unsigned int r = 0; r < 1; r++)
    {
      rr = (r + 1)*(r + 1 + 1)*(r + 1 + 2)/6 + 1*(1 + 1)/2;
      ss = r*(r + 1)*(r + 2)/6;
      basisvalues[rr] = basisvalues[ss]*(r*(1.000000000000000 + Y) + (2.000000000000000 + Z + 3.000000000000000*Y)/2.000000000000000);
    }// end loop over 'r'
    for (unsigned int r = 0; r < 1; r++)
    {
      for (unsigned int s = 0; s < 1 - r; s++)
      {
        rr = (r + s + 1)*(r + s + 1 + 1)*(r + s + 1 + 2)/6 + (s + 1)*(s + 1 + 1)/2 + 1;
        ss = (r + s)*(r + s + 1)*(r + s + 2)/6 + s*(s + 1)/2;
        basisvalues[rr] = basisvalues[ss]*(1.000000000000000 + r + s + Z*(2.000000000000000 + r + s));
      }// end loop over 's'
    }// end loop over 'r'
    for (unsigned int r = 0; r < 2; r++)
    {
      for (unsigned int s = 0; s < 2 - r; s++)
      {
        for (unsigned int t = 0; t < 2 - r - s; t++)
        {
          rr = (r + s + t)*(r + s + t + 1)*(r + s + t + 2)/6 + (s + t)*(s + t + 1)/2 + t;
          basisvalues[rr] *= std::sqrt((0.500000000000000 + r)*(1.000000000000000 + r + s)*(1.500000000000000 + r + s + t));
        }// end loop over 't'
      }// end loop over 's'
    }// end loop over 'r'
    
    // Table(s) of coefficients.
    static const double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.000000000000000, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0.000000000000000, 0.000000000000000, 0.223606797749979}};
    
    // Compute value(s).
    for (unsigned int r = 0; r < 4; r++)
    {
      *values += coefficients0[dof][r]*basisvalues[r];
    }// end loop over 'r'
  }

  /// Evaluate all basis functions at given point in cell
  virtual void evaluate_basis_all(double* values,
                                  const double* coordinates,
                                  const ufc::cell& c) const
  {
    // Helper variable to hold values of a single dof.
    double dof_values = 0.000000000000000;
    
    // Loop dofs and call evaluate_basis.
    for (unsigned int r = 0; r < 4; r++)
    {
      evaluate_basis(r, &dof_values, coordinates, c);
      values[r] = dof_values;
    }// end loop over 'r'
  }

  /// Evaluate order n derivatives of basis function i at given point in cell
  virtual void evaluate_basis_derivatives(unsigned int i,
                                          unsigned int n,
                                          double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const
  {
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_02 = x[3][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    const double J_12 = x[3][1] - x[0][1];
    const double J_20 = x[1][2] - x[0][2];
    const double J_21 = x[2][2] - x[0][2];
    const double J_22 = x[3][2] - x[0][2];
    
    // Compute sub determinants
    const double d_00 = J_11*J_22 - J_12*J_21;
    const double d_01 = J_12*J_20 - J_10*J_22;
    const double d_02 = J_10*J_21 - J_11*J_20;
    const double d_10 = J_02*J_21 - J_01*J_22;
    const double d_11 = J_00*J_22 - J_02*J_20;
    const double d_12 = J_01*J_20 - J_00*J_21;
    const double d_20 = J_01*J_12 - J_02*J_11;
    const double d_21 = J_02*J_10 - J_00*J_12;
    const double d_22 = J_00*J_11 - J_01*J_10;
    
    // Compute determinant of Jacobian
    double detJ = J_00*d_00 + J_10*d_10 + J_20*d_20;
    
    // Compute inverse of Jacobian
    const double K_00 = d_00 / detJ;
    const double K_01 = d_10 / detJ;
    const double K_02 = d_20 / detJ;
    const double K_10 = d_01 / detJ;
    const double K_11 = d_11 / detJ;
    const double K_12 = d_21 / detJ;
    const double K_20 = d_02 / detJ;
    const double K_21 = d_12 / detJ;
    const double K_22 = d_22 / detJ;
    
    // Compute constants
    const double C0 = x[3][0] + x[2][0] + x[1][0] - x[0][0];
    const double C1 = x[3][1] + x[2][1] + x[1][1] - x[0][1];
    const double C2 = x[3][2] + x[2][2] + x[1][2] - x[0][2];
    
    // Get coordinates and map to the reference (FIAT) element
    double X = (d_00*(2.0*coordinates[0] - C0) + d_10*(2.0*coordinates[1] - C1) + d_20*(2.0*coordinates[2] - C2)) / detJ;
    double Y = (d_01*(2.0*coordinates[0] - C0) + d_11*(2.0*coordinates[1] - C1) + d_21*(2.0*coordinates[2] - C2)) / detJ;
    double Z = (d_02*(2.0*coordinates[0] - C0) + d_12*(2.0*coordinates[1] - C1) + d_22*(2.0*coordinates[2] - C2)) / detJ;
    
    
    // Compute number of derivatives.
    unsigned int num_derivatives = 1;
    for (unsigned int r = 0; r < n; r++)
    {
      num_derivatives *= 3;
    }// end loop over 'r'
    
    // Declare pointer to two dimensional array that holds combinations of derivatives and initialise
    unsigned int **combinations = new unsigned int *[num_derivatives];
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      combinations[row] = new unsigned int [n];
      for (unsigned int col = 0; col < n; col++)
        combinations[row][col] = 0;
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
    const double Jinv[3][3] = {{K_00, K_01, K_02}, {K_10, K_11, K_12}, {K_20, K_21, K_22}};
    
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
    
    // Reset values. Assuming that values is always an array.
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      values[r] = 0.000000000000000;
    }// end loop over 'r'
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Array of basisvalues.
    double basisvalues[4] = {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000};
    
    // Declare helper variables.
    unsigned int rr = 0;
    unsigned int ss = 0;
    double tmp0 = 0.500000000000000*(2.000000000000000 + Y + Z + 2.000000000000000*X);
    
    // Compute basisvalues.
    basisvalues[0] = 1.000000000000000;
    basisvalues[1] = tmp0;
    for (unsigned int r = 0; r < 1; r++)
    {
      rr = (r + 1)*(r + 1 + 1)*(r + 1 + 2)/6 + 1*(1 + 1)/2;
      ss = r*(r + 1)*(r + 2)/6;
      basisvalues[rr] = basisvalues[ss]*(r*(1.000000000000000 + Y) + (2.000000000000000 + Z + 3.000000000000000*Y)/2.000000000000000);
    }// end loop over 'r'
    for (unsigned int r = 0; r < 1; r++)
    {
      for (unsigned int s = 0; s < 1 - r; s++)
      {
        rr = (r + s + 1)*(r + s + 1 + 1)*(r + s + 1 + 2)/6 + (s + 1)*(s + 1 + 1)/2 + 1;
        ss = (r + s)*(r + s + 1)*(r + s + 2)/6 + s*(s + 1)/2;
        basisvalues[rr] = basisvalues[ss]*(1.000000000000000 + r + s + Z*(2.000000000000000 + r + s));
      }// end loop over 's'
    }// end loop over 'r'
    for (unsigned int r = 0; r < 2; r++)
    {
      for (unsigned int s = 0; s < 2 - r; s++)
      {
        for (unsigned int t = 0; t < 2 - r - s; t++)
        {
          rr = (r + s + t)*(r + s + t + 1)*(r + s + t + 2)/6 + (s + t)*(s + t + 1)/2 + t;
          basisvalues[rr] *= std::sqrt((0.500000000000000 + r)*(1.000000000000000 + r + s)*(1.500000000000000 + r + s + t));
        }// end loop over 't'
      }// end loop over 's'
    }// end loop over 'r'
    
    // Table(s) of coefficients.
    static const double coefficients0[4][4] = \
    {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
    {0.288675134594813, 0.000000000000000, 0.210818510677892, -0.074535599249993},
    {0.288675134594813, 0.000000000000000, 0.000000000000000, 0.223606797749979}};
    
    // Tables of derivatives of the polynomial base (transpose).
    static const double dmats0[4][4] = \
    {{0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {6.324555320336760, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000}};
    
    static const double dmats1[4][4] = \
    {{0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {3.162277660168380, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {5.477225575051662, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000}};
    
    static const double dmats2[4][4] = \
    {{0.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {3.162277660168380, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {1.825741858350554, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {5.163977794943223, 0.000000000000000, 0.000000000000000, 0.000000000000000}};
    
    // Compute reference derivatives.
    // Declare pointer to array of derivatives on FIAT element.
    double *derivatives = new double[num_derivatives];
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      derivatives[r] = 0.000000000000000;
    }// end loop over 'r'
    
    // Declare derivative matrix (of polynomial basis).
    double dmats[4][4] = \
    {{1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000}};
    
    // Declare (auxiliary) derivative matrix (of polynomial basis).
    double dmats_old[4][4] = \
    {{1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000},
    {0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000}};
    
    // Loop possible derivatives.
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      // Resetting dmats values to compute next derivative.
      for (unsigned int t = 0; t < 4; t++)
      {
        for (unsigned int u = 0; u < 4; u++)
        {
          dmats[t][u] = 0.000000000000000;
          if (t == u)
          {
          dmats[t][u] = 1.000000000000000;
          }
          
        }// end loop over 'u'
      }// end loop over 't'
      
      // Looping derivative order to generate dmats.
      for (unsigned int s = 0; s < n; s++)
      {
        // Updating dmats_old with new values and resetting dmats.
        for (unsigned int t = 0; t < 4; t++)
        {
          for (unsigned int u = 0; u < 4; u++)
          {
            dmats_old[t][u] = dmats[t][u];
            dmats[t][u] = 0.000000000000000;
          }// end loop over 'u'
        }// end loop over 't'
        
        // Update dmats using an inner product.
        if (combinations[r][s] == 0)
        {
        for (unsigned int t = 0; t < 4; t++)
        {
          for (unsigned int u = 0; u < 4; u++)
          {
            for (unsigned int tu = 0; tu < 4; tu++)
            {
              dmats[t][u] += dmats0[t][tu]*dmats_old[tu][u];
            }// end loop over 'tu'
          }// end loop over 'u'
        }// end loop over 't'
        }
        
        if (combinations[r][s] == 1)
        {
        for (unsigned int t = 0; t < 4; t++)
        {
          for (unsigned int u = 0; u < 4; u++)
          {
            for (unsigned int tu = 0; tu < 4; tu++)
            {
              dmats[t][u] += dmats1[t][tu]*dmats_old[tu][u];
            }// end loop over 'tu'
          }// end loop over 'u'
        }// end loop over 't'
        }
        
        if (combinations[r][s] == 2)
        {
        for (unsigned int t = 0; t < 4; t++)
        {
          for (unsigned int u = 0; u < 4; u++)
          {
            for (unsigned int tu = 0; tu < 4; tu++)
            {
              dmats[t][u] += dmats2[t][tu]*dmats_old[tu][u];
            }// end loop over 'tu'
          }// end loop over 'u'
        }// end loop over 't'
        }
        
      }// end loop over 's'
      for (unsigned int s = 0; s < 4; s++)
      {
        for (unsigned int t = 0; t < 4; t++)
        {
          derivatives[r] += coefficients0[dof][s]*dmats[s][t]*basisvalues[t];
        }// end loop over 't'
      }// end loop over 's'
    }// end loop over 'r'
    
    // Transform derivatives back to physical element
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      for (unsigned int s = 0; s < num_derivatives; s++)
      {
        values[r] += transform[r][s]*derivatives[s];
      }// end loop over 's'
    }// end loop over 'r'
    
    // Delete pointer to array of derivatives on FIAT element
    delete [] derivatives;
    
    // Delete pointer to array of combinations of derivatives and transform
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      delete [] combinations[r];
    }// end loop over 'r'
    delete [] combinations;
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      delete [] transform[r];
    }// end loop over 'r'
    delete [] transform;
  }

  /// Evaluate order n derivatives of all basis functions at given point in cell
  virtual void evaluate_basis_derivatives_all(unsigned int n,
                                              double* values,
                                              const double* coordinates,
                                              const ufc::cell& c) const
  {
    // Compute number of derivatives.
    unsigned int num_derivatives = 1;
    for (unsigned int r = 0; r < n; r++)
    {
      num_derivatives *= 3;
    }// end loop over 'r'
    
    // Helper variable to hold values of a single dof.
    double *dof_values = new double[num_derivatives];
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      dof_values[r] = 0.000000000000000;
    }// end loop over 'r'
    
    // Loop dofs and call evaluate_basis_derivatives.
    for (unsigned int r = 0; r < 4; r++)
    {
      evaluate_basis_derivatives(r, n, dof_values, coordinates, c);
      for (unsigned int s = 0; s < num_derivatives; s++)
      {
        values[r*num_derivatives + s] = dof_values[s];
      }// end loop over 's'
    }// end loop over 'r'
    
    // Delete pointer.
    delete [] dof_values;
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    // Declare variables for result of evaluation.
    double vals[1];
    
    // Declare variable for physical coordinates.
    double y[3];
    const double * const * x = c.coordinates;
    switch (i)
    {
    case 0:
      {
        y[0] = x[0][0];
      y[1] = x[0][1];
      y[2] = x[0][2];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    case 1:
      {
        y[0] = x[1][0];
      y[1] = x[1][1];
      y[2] = x[1][2];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    case 2:
      {
        y[0] = x[2][0];
      y[1] = x[2][1];
      y[2] = x[2][2];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    case 3:
      {
        y[0] = x[3][0];
      y[1] = x[3][1];
      y[2] = x[3][2];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    }
    
    return 0.000000000000000;
  }

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const ufc::cell& c) const
  {
    // Declare variables for result of evaluation.
    double vals[1];
    
    // Declare variable for physical coordinates.
    double y[3];
    const double * const * x = c.coordinates;
    y[0] = x[0][0];
    y[1] = x[0][1];
    y[2] = x[0][2];
    f.evaluate(vals, y, c);
    values[0] = vals[0];
    y[0] = x[1][0];
    y[1] = x[1][1];
    y[2] = x[1][2];
    f.evaluate(vals, y, c);
    values[1] = vals[0];
    y[0] = x[2][0];
    y[1] = x[2][1];
    y[2] = x[2][2];
    f.evaluate(vals, y, c);
    values[2] = vals[0];
    y[0] = x[3][0];
    y[1] = x[3][1];
    y[2] = x[3][2];
    f.evaluate(vals, y, c);
    values[3] = vals[0];
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
    // Evaluate function and change variables
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
    vertex_values[3] = dof_values[3];
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 0;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return 0;
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class massmatrix3d_dof_map_0: public ufc::dof_map
{
private:

  unsigned int _global_dimension;
public:

  /// Constructor
  massmatrix3d_dof_map_0() : ufc::dof_map()
  {
    _global_dimension = 0;
  }

  /// Destructor
  virtual ~massmatrix3d_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dofmap for FiniteElement('Lagrange', Cell('tetrahedron', 1, Space(3)), 1)";
  }

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const
  {
    switch (d)
    {
    case 0:
      {
        return true;
        break;
      }
    case 1:
      {
        return false;
        break;
      }
    case 2:
      {
        return false;
        break;
      }
    case 3:
      {
        return false;
        break;
      }
    }
    
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    _global_dimension = m.num_entities[0];
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
    return _global_dimension;
  }

  /// Return the dimension of the local finite element function space for a cell
  virtual unsigned int local_dimension(const ufc::cell& c) const
  {
    return 4;
  }

  /// Return the maximum dimension of the local finite element function space
  virtual unsigned int max_local_dimension() const
  {
    return 4;
  }

  // Return the geometric dimension of the coordinates this dof map provides
  virtual unsigned int geometric_dimension() const
  {
    return 3;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 3;
  }

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual unsigned int num_entity_dofs(unsigned int d) const
  {
    switch (d)
    {
    case 0:
      {
        return 1;
        break;
      }
    case 1:
      {
        return 0;
        break;
      }
    case 2:
      {
        return 0;
        break;
      }
    case 3:
      {
        return 0;
        break;
      }
    }
    
    return 0;
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
                                   unsigned int facet) const
  {
    switch (facet)
    {
    case 0:
      {
        dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
        break;
      }
    case 1:
      {
        dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
        break;
      }
    case 2:
      {
        dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
        break;
      }
    case 3:
      {
        dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
        break;
      }
    }
    
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(unsigned int* dofs,
                                    unsigned int d, unsigned int i) const
  {
    if (d > 3)
    {
    throw std::runtime_error("d is larger than dimension (3)");
    }
    
    switch (d)
    {
    case 0:
      {
        if (i > 3)
      {
      throw std::runtime_error("i is larger than number of entities (3)");
      }
      
      switch (i)
      {
      case 0:
        {
          dofs[0] = 0;
          break;
        }
      case 1:
        {
          dofs[0] = 1;
          break;
        }
      case 2:
        {
          dofs[0] = 2;
          break;
        }
      case 3:
        {
          dofs[0] = 3;
          break;
        }
      }
      
        break;
      }
    case 1:
      {
        
        break;
      }
    case 2:
      {
        
        break;
      }
    case 3:
      {
        
        break;
      }
    }
    
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
    return 0;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return 0;
  }

};

/// This class defines the interface for the tabulation of the cell
/// tensor corresponding to the local contribution to a form from
/// the integral over a cell.

class massmatrix3d_cell_integral_0_0: public ufc::cell_integral
{
public:

  /// Constructor
  massmatrix3d_cell_integral_0_0() : ufc::cell_integral()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~massmatrix3d_cell_integral_0_0()
  {
    // Do nothing
  }

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c) const
  {
    // Number of operations (multiply-add pairs) for Jacobian data:      18
    // Number of operations (multiply-add pairs) for geometry tensor:    0
    // Number of operations (multiply-add pairs) for tensor contraction: 8
    // Total number of operations (multiply-add pairs):                  26
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    const double J_01 = x[2][0] - x[0][0];
    const double J_02 = x[3][0] - x[0][0];
    const double J_10 = x[1][1] - x[0][1];
    const double J_11 = x[2][1] - x[0][1];
    const double J_12 = x[3][1] - x[0][1];
    const double J_20 = x[1][2] - x[0][2];
    const double J_21 = x[2][2] - x[0][2];
    const double J_22 = x[3][2] - x[0][2];
    
    // Compute sub determinants
    const double d_00 = J_11*J_22 - J_12*J_21;
    const double d_10 = J_02*J_21 - J_01*J_22;
    const double d_20 = J_01*J_12 - J_02*J_11;
    
    // Compute determinant of Jacobian
    double detJ = J_00*d_00 + J_10*d_10 + J_20*d_20;
    
    // Compute inverse of Jacobian
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    // Compute geometry tensor
    const double G0_ = det;
    
    // Compute element tensor
    A[0] = 0.016666666666667*G0_;
    A[1] = 0.008333333333333*G0_;
    A[2] = 0.008333333333333*G0_;
    A[3] = 0.008333333333333*G0_;
    A[4] = 0.008333333333333*G0_;
    A[5] = 0.016666666666667*G0_;
    A[6] = 0.008333333333333*G0_;
    A[7] = 0.008333333333333*G0_;
    A[8] = 0.008333333333333*G0_;
    A[9] = 0.008333333333333*G0_;
    A[10] = 0.016666666666667*G0_;
    A[11] = 0.008333333333333*G0_;
    A[12] = 0.008333333333333*G0_;
    A[13] = 0.008333333333333*G0_;
    A[14] = 0.008333333333333*G0_;
    A[15] = 0.016666666666667*G0_;
  }

};

/// This class defines the interface for the assembly of the global
/// tensor corresponding to a form with r + n arguments, that is, a
/// mapping
///
///     a : V1 x V2 x ... Vr x W1 x W2 x ... x Wn -> R
///
/// with arguments v1, v2, ..., vr, w1, w2, ..., wn. The rank r
/// global tensor A is defined by
///
///     A = a(V1, V2, ..., Vr, w1, w2, ..., wn),
///
/// where each argument Vj represents the application to the
/// sequence of basis functions of Vj and w1, w2, ..., wn are given
/// fixed functions (coefficients).

class massmatrix3d_form_0: public ufc::form
{
public:

  /// Constructor
  massmatrix3d_form_0() : ufc::form()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~massmatrix3d_form_0()
  {
    // Do nothing
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
    return "Form([Integral(Product(Argument(FiniteElement('Lagrange', Cell('tetrahedron', 1, Space(3)), 1), 0), Argument(FiniteElement('Lagrange', Cell('tetrahedron', 1, Space(3)), 1), 1)), Measure('cell', 0, None))])";
  }

  /// Return the rank of the global tensor (r)
  virtual unsigned int rank() const
  {
    return 2;
  }

  /// Return the number of coefficients (n)
  virtual unsigned int num_coefficients() const
  {
    return 0;
  }

  /// Return the number of cell integrals
  virtual unsigned int num_cell_integrals() const
  {
    return 1;
  }

  /// Return the number of exterior facet integrals
  virtual unsigned int num_exterior_facet_integrals() const
  {
    return 0;
  }

  /// Return the number of interior facet integrals
  virtual unsigned int num_interior_facet_integrals() const
  {
    return 0;
  }

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(unsigned int i) const
  {
    switch (i)
    {
    case 0:
      {
        return new massmatrix3d_finite_element_0();
        break;
      }
    case 1:
      {
        return new massmatrix3d_finite_element_0();
        break;
      }
    }
    
    return 0;
  }

  /// Create a new dof map for argument function i
  virtual ufc::dof_map* create_dof_map(unsigned int i) const
  {
    switch (i)
    {
    case 0:
      {
        return new massmatrix3d_dof_map_0();
        break;
      }
    case 1:
      {
        return new massmatrix3d_dof_map_0();
        break;
      }
    }
    
    return 0;
  }

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(unsigned int i) const
  {
    switch (i)
    {
    case 0:
      {
        return new massmatrix3d_cell_integral_0_0();
        break;
      }
    }
    
    return 0;
  }

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(unsigned int i) const
  {
    return 0;
  }

  /// Create a new interior facet integral on sub domain i
  virtual ufc::interior_facet_integral* create_interior_facet_integral(unsigned int i) const
  {
    return 0;
  }

};

// DOLFIN wrappers

// Standard library includes
#include <string>

// DOLFIN includes
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/Form.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/CoefficientAssigner.h>

namespace MassMatrix3D
{

class Form_0_FunctionSpace_0: public dolfin::FunctionSpace
{
public:

  Form_0_FunctionSpace_0(const dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_0(dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_0(boost::shared_ptr<dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), *mesh)))
  {
      // Do nothing
  }

  Form_0_FunctionSpace_0(boost::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), *mesh)))
  {
      // Do nothing
  }

  ~Form_0_FunctionSpace_0()
  {
  }

};

class Form_0_FunctionSpace_1: public dolfin::FunctionSpace
{
public:

  Form_0_FunctionSpace_1(const dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_1(dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_1(boost::shared_ptr<dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), *mesh)))
  {
      // Do nothing
  }

  Form_0_FunctionSpace_1(boost::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new massmatrix3d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new massmatrix3d_dof_map_0()), *mesh)))
  {
      // Do nothing
  }

  ~Form_0_FunctionSpace_1()
  {
  }

};

class Form_0: public dolfin::Form
{
public:

  // Constructor
  Form_0(const dolfin::FunctionSpace& V0, const dolfin::FunctionSpace& V1):
    dolfin::Form(2, 0)
  {
    _function_spaces[0] = reference_to_no_delete_pointer(V0);
    _function_spaces[1] = reference_to_no_delete_pointer(V1);

    _ufc_form = boost::shared_ptr<const ufc::form>(new massmatrix3d_form_0());
  }

  // Constructor
  Form_0(boost::shared_ptr<const dolfin::FunctionSpace> V0, boost::shared_ptr<const dolfin::FunctionSpace> V1):
    dolfin::Form(2, 0)
  {
    _function_spaces[0] = V0;
    _function_spaces[1] = V1;

    _ufc_form = boost::shared_ptr<const ufc::form>(new massmatrix3d_form_0());
  }

  // Destructor
  ~Form_0()
  {}

  /// Return the number of the coefficient with this name
  virtual dolfin::uint coefficient_number(const std::string& name) const
  {

    dolfin::error("No coefficients.");
    return 0;
  }

  /// Return the name of the coefficient with this number
  virtual std::string coefficient_name(dolfin::uint i) const
  {

    dolfin::error("No coefficients.");
    return "unnamed";
  }

  // Typedefs
  typedef Form_0_FunctionSpace_0 TestSpace;
  typedef Form_0_FunctionSpace_1 TrialSpace;

  // Coefficients
};

// Class typedefs
typedef Form_0 BilinearForm;
typedef Form_0::TestSpace FunctionSpace;

}

#endif
