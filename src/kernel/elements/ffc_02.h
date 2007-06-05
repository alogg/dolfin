// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.4.0-pre1.

#ifndef __FFC_02_H
#define __FFC_02_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_02_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_02_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_02_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Lagrange finite element of degree 3 on a triangle";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::triangle;
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
    return 10;
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
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
      
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute constants
    const double C0 = element_coordinates[1][0] + element_coordinates[2][0];
    const double C1 = element_coordinates[1][1] + element_coordinates[2][1];
    
    // Get coordinates and map to the reference (FIAT) element
    double x = (J_01*C1 - J_11*C0 + 2.0*J_11*coordinates[0] - 2.0*J_01*coordinates[1]) / detJ;
    double y = (J_10*C0 - J_00*C1 - 2.0*J_10*coordinates[0] + 2.0*J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 * (1.0 + x)/(1.0 - y) - 1.0;
    
    // Reset values
    *values = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    const double scalings_y_3 = scalings_y_2*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    const double psitilde_a_3 = 1.66666666666667*x*psitilde_a_2 - 0.666666666666667*psitilde_a_1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_0_3 = 0.05*psitilde_bs_0_2 + 1.75*y*psitilde_bs_0_2 - 0.7*psitilde_bs_0_1;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_1_2 = 0.54*psitilde_bs_1_1 + 2.1*y*psitilde_bs_1_1 - 0.56*psitilde_bs_1_0;
    const double psitilde_bs_2_0 = 1;
    const double psitilde_bs_2_1 = 3.5*y + 2.5;
    const double psitilde_bs_3_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    const double basisvalue6 = 3.74165738677394*psitilde_a_3*scalings_y_3*psitilde_bs_3_0;
    const double basisvalue7 = 3.16227766016838*psitilde_a_2*scalings_y_2*psitilde_bs_2_1;
    const double basisvalue8 = 2.44948974278318*psitilde_a_1*scalings_y_1*psitilde_bs_1_2;
    const double basisvalue9 = 1.4142135623731*psitilde_a_0*scalings_y_0*psitilde_bs_0_3;
    
    // Table(s) of coefficients
    const static double coefficients0[10][10] = \
    {{0.0471404520791032, -0.0288675134594813, -0.0166666666666667, 0.0782460796435952, 0.0606091526731326, 0.0349927106111883, -0.0601337794302955, -0.0508223195384204, -0.0393667994375868, -0.0227284322524247},
    {0.0471404520791031, 0.0288675134594813, -0.0166666666666667, 0.0782460796435952, -0.0606091526731326, 0.0349927106111883, 0.0601337794302955, -0.0508223195384204, 0.0393667994375868, -0.0227284322524248},
    {0.0471404520791032, 0, 0.0333333333333333, 0, 0, 0.104978131833565, 0, 0, 0, 0.090913729009699},
    {0.106066017177982, 0.259807621135332, -0.15, 0.117369119465393, 0.0606091526731326, -0.0787335988751736, 0, 0.101644639076841, -0.131222664791956, 0.090913729009699},
    {0.106066017177982, 0, 0.3, 0, 0.151522881682832, 0.0262445329583912, 0, 0, 0.131222664791956, -0.136370593514548},
    {0.106066017177982, 0, 0.3, 0, -0.151522881682832, 0.0262445329583912, 0, 0, -0.131222664791956, -0.136370593514548},
    {0.106066017177982, -0.259807621135332, -0.15, 0.117369119465393, -0.0606091526731327, -0.0787335988751736, 0, 0.101644639076841, 0.131222664791956, 0.0909137290096989},
    {0.106066017177982, -0.259807621135332, -0.15, -0.0782460796435952, 0.090913729009699, 0.0962299541807677, 0.180401338290886, 0.0508223195384204, -0.0131222664791956, -0.0227284322524248},
    {0.106066017177982, 0.259807621135332, -0.15, -0.0782460796435952, -0.090913729009699, 0.0962299541807677, -0.180401338290886, 0.0508223195384204, 0.0131222664791956, -0.0227284322524247},
    {0.636396103067893, 0, 0, -0.234738238930785, 0, -0.262445329583912, 0, -0.203289278153682, 0, 0.090913729009699}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    const double coeff0_4 = coefficients0[dof][4];
    const double coeff0_5 = coefficients0[dof][5];
    const double coeff0_6 = coefficients0[dof][6];
    const double coeff0_7 = coefficients0[dof][7];
    const double coeff0_8 = coefficients0[dof][8];
    const double coeff0_9 = coefficients0[dof][9];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5 + coeff0_6*basisvalue6 + coeff0_7*basisvalue7 + coeff0_8*basisvalue8 + coeff0_9*basisvalue9;
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
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
      
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute constants
    const double C0 = element_coordinates[1][0] + element_coordinates[2][0];
    const double C1 = element_coordinates[1][1] + element_coordinates[2][1];
    
    // Get coordinates and map to the reference (FIAT) element
    double x = (J_01*C1 - J_11*C0 + 2.0*J_11*coordinates[0] - 2.0*J_01*coordinates[1]) / detJ;
    double y = (J_10*C0 - J_00*C1 - 2.0*J_10*coordinates[0] + 2.0*J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 * (1.0 + x)/(1.0 - y) - 1.0;
    
    // Compute number of derivatives
    unsigned int num_derivatives = 1;
    
    for (unsigned int j = 0; j < n; j++)
      num_derivatives *= 2;
    
    
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
          if (combinations[row][col] + 1 > 1)
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
    const double Jinv[2][2] =  {{2*J_11 / detJ, -2*J_01 / detJ}, {-2*J_10 / detJ, 2*J_00 / detJ}};
    
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
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    const double scalings_y_3 = scalings_y_2*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    const double psitilde_a_3 = 1.66666666666667*x*psitilde_a_2 - 0.666666666666667*psitilde_a_1;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_0_3 = 0.05*psitilde_bs_0_2 + 1.75*y*psitilde_bs_0_2 - 0.7*psitilde_bs_0_1;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_1_2 = 0.54*psitilde_bs_1_1 + 2.1*y*psitilde_bs_1_1 - 0.56*psitilde_bs_1_0;
    const double psitilde_bs_2_0 = 1;
    const double psitilde_bs_2_1 = 3.5*y + 2.5;
    const double psitilde_bs_3_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    const double basisvalue6 = 3.74165738677394*psitilde_a_3*scalings_y_3*psitilde_bs_3_0;
    const double basisvalue7 = 3.16227766016838*psitilde_a_2*scalings_y_2*psitilde_bs_2_1;
    const double basisvalue8 = 2.44948974278318*psitilde_a_1*scalings_y_1*psitilde_bs_1_2;
    const double basisvalue9 = 1.4142135623731*psitilde_a_0*scalings_y_0*psitilde_bs_0_3;
    
    // Table(s) of coefficients
    const static double coefficients0[10][10] = \
    {{0.0471404520791032, -0.0288675134594813, -0.0166666666666667, 0.0782460796435952, 0.0606091526731326, 0.0349927106111883, -0.0601337794302955, -0.0508223195384204, -0.0393667994375868, -0.0227284322524247},
    {0.0471404520791031, 0.0288675134594813, -0.0166666666666667, 0.0782460796435952, -0.0606091526731326, 0.0349927106111883, 0.0601337794302955, -0.0508223195384204, 0.0393667994375868, -0.0227284322524248},
    {0.0471404520791032, 0, 0.0333333333333333, 0, 0, 0.104978131833565, 0, 0, 0, 0.090913729009699},
    {0.106066017177982, 0.259807621135332, -0.15, 0.117369119465393, 0.0606091526731326, -0.0787335988751736, 0, 0.101644639076841, -0.131222664791956, 0.090913729009699},
    {0.106066017177982, 0, 0.3, 0, 0.151522881682832, 0.0262445329583912, 0, 0, 0.131222664791956, -0.136370593514548},
    {0.106066017177982, 0, 0.3, 0, -0.151522881682832, 0.0262445329583912, 0, 0, -0.131222664791956, -0.136370593514548},
    {0.106066017177982, -0.259807621135332, -0.15, 0.117369119465393, -0.0606091526731327, -0.0787335988751736, 0, 0.101644639076841, 0.131222664791956, 0.0909137290096989},
    {0.106066017177982, -0.259807621135332, -0.15, -0.0782460796435952, 0.090913729009699, 0.0962299541807677, 0.180401338290886, 0.0508223195384204, -0.0131222664791956, -0.0227284322524248},
    {0.106066017177982, 0.259807621135332, -0.15, -0.0782460796435952, -0.090913729009699, 0.0962299541807677, -0.180401338290886, 0.0508223195384204, 0.0131222664791956, -0.0227284322524247},
    {0.636396103067893, 0, 0, -0.234738238930785, 0, -0.262445329583912, 0, -0.203289278153682, 0, 0.090913729009699}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.44948974278318, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 4.74341649025257, 0, 0, 0, 0, 0, 0, 0, 0},
    {2, 0, 3.53553390593274, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.64575131106459, 0, -1.49666295470958, 6.83130051063973, 0, 0.30550504633039, 0, 0, 0, 0},
    {0, 2.19089023002067, 0, 0, 6.26099033699941, 0, 0, 0, 0, 0},
    {1.73205080756888, 0, 3.91918358845308, 0, 0, 4.2, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    const static double dmats1[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.22474487139159, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.12132034355964, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.29099444873581, 2.37170824512628, -0.456435464587638, 0, 0, 0, 0, 0, 0, 0},
    {1, 3.06186217847897, 1.76776695296637, 0, 0, 0, 0, 0, 0, 0},
    {-1.15470053837925, 0, 4.08248290463863, 0, 0, 0, 0, 0, 0, 0},
    {1.3228756555323, 2.59229627936314, -0.748331477354788, 3.41565025531987, -0.529150262212918, 0.152752523165196, 0, 0, 0, 0},
    {1.11803398874989, 1.09544511501033, 1.26491106406735, 4.04145188432738, 3.13049516849971, -0.903696114115064, 0, 0, 0, 0},
    {0.866025403784438, -2.54558441227157, 1.95959179422654, 0, 4.84974226119286, 2.1, 0, 0, 0, 0},
    {2.5, 0, -1.4142135623731, 0, 0, 6.06217782649107, 0, 0, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff0_3 = 0;
    double coeff0_4 = 0;
    double coeff0_5 = 0;
    double coeff0_6 = 0;
    double coeff0_7 = 0;
    double coeff0_8 = 0;
    double coeff0_9 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    double new_coeff0_4 = 0;
    double new_coeff0_5 = 0;
    double new_coeff0_6 = 0;
    double new_coeff0_7 = 0;
    double new_coeff0_8 = 0;
    double new_coeff0_9 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff0_3 = coefficients0[dof][3];
      new_coeff0_4 = coefficients0[dof][4];
      new_coeff0_5 = coefficients0[dof][5];
      new_coeff0_6 = coefficients0[dof][6];
      new_coeff0_7 = coefficients0[dof][7];
      new_coeff0_8 = coefficients0[dof][8];
      new_coeff0_9 = coefficients0[dof][9];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff0_3 = new_coeff0_3;
        coeff0_4 = new_coeff0_4;
        coeff0_5 = new_coeff0_5;
        coeff0_6 = new_coeff0_6;
        coeff0_7 = new_coeff0_7;
        coeff0_8 = new_coeff0_8;
        coeff0_9 = new_coeff0_9;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0] + coeff0_6*dmats0[6][0] + coeff0_7*dmats0[7][0] + coeff0_8*dmats0[8][0] + coeff0_9*dmats0[9][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1] + coeff0_6*dmats0[6][1] + coeff0_7*dmats0[7][1] + coeff0_8*dmats0[8][1] + coeff0_9*dmats0[9][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2] + coeff0_6*dmats0[6][2] + coeff0_7*dmats0[7][2] + coeff0_8*dmats0[8][2] + coeff0_9*dmats0[9][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3] + coeff0_6*dmats0[6][3] + coeff0_7*dmats0[7][3] + coeff0_8*dmats0[8][3] + coeff0_9*dmats0[9][3];
          new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4] + coeff0_6*dmats0[6][4] + coeff0_7*dmats0[7][4] + coeff0_8*dmats0[8][4] + coeff0_9*dmats0[9][4];
          new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5] + coeff0_6*dmats0[6][5] + coeff0_7*dmats0[7][5] + coeff0_8*dmats0[8][5] + coeff0_9*dmats0[9][5];
          new_coeff0_6 = coeff0_0*dmats0[0][6] + coeff0_1*dmats0[1][6] + coeff0_2*dmats0[2][6] + coeff0_3*dmats0[3][6] + coeff0_4*dmats0[4][6] + coeff0_5*dmats0[5][6] + coeff0_6*dmats0[6][6] + coeff0_7*dmats0[7][6] + coeff0_8*dmats0[8][6] + coeff0_9*dmats0[9][6];
          new_coeff0_7 = coeff0_0*dmats0[0][7] + coeff0_1*dmats0[1][7] + coeff0_2*dmats0[2][7] + coeff0_3*dmats0[3][7] + coeff0_4*dmats0[4][7] + coeff0_5*dmats0[5][7] + coeff0_6*dmats0[6][7] + coeff0_7*dmats0[7][7] + coeff0_8*dmats0[8][7] + coeff0_9*dmats0[9][7];
          new_coeff0_8 = coeff0_0*dmats0[0][8] + coeff0_1*dmats0[1][8] + coeff0_2*dmats0[2][8] + coeff0_3*dmats0[3][8] + coeff0_4*dmats0[4][8] + coeff0_5*dmats0[5][8] + coeff0_6*dmats0[6][8] + coeff0_7*dmats0[7][8] + coeff0_8*dmats0[8][8] + coeff0_9*dmats0[9][8];
          new_coeff0_9 = coeff0_0*dmats0[0][9] + coeff0_1*dmats0[1][9] + coeff0_2*dmats0[2][9] + coeff0_3*dmats0[3][9] + coeff0_4*dmats0[4][9] + coeff0_5*dmats0[5][9] + coeff0_6*dmats0[6][9] + coeff0_7*dmats0[7][9] + coeff0_8*dmats0[8][9] + coeff0_9*dmats0[9][9];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0] + coeff0_6*dmats1[6][0] + coeff0_7*dmats1[7][0] + coeff0_8*dmats1[8][0] + coeff0_9*dmats1[9][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1] + coeff0_6*dmats1[6][1] + coeff0_7*dmats1[7][1] + coeff0_8*dmats1[8][1] + coeff0_9*dmats1[9][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2] + coeff0_6*dmats1[6][2] + coeff0_7*dmats1[7][2] + coeff0_8*dmats1[8][2] + coeff0_9*dmats1[9][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3] + coeff0_6*dmats1[6][3] + coeff0_7*dmats1[7][3] + coeff0_8*dmats1[8][3] + coeff0_9*dmats1[9][3];
          new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4] + coeff0_6*dmats1[6][4] + coeff0_7*dmats1[7][4] + coeff0_8*dmats1[8][4] + coeff0_9*dmats1[9][4];
          new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5] + coeff0_6*dmats1[6][5] + coeff0_7*dmats1[7][5] + coeff0_8*dmats1[8][5] + coeff0_9*dmats1[9][5];
          new_coeff0_6 = coeff0_0*dmats1[0][6] + coeff0_1*dmats1[1][6] + coeff0_2*dmats1[2][6] + coeff0_3*dmats1[3][6] + coeff0_4*dmats1[4][6] + coeff0_5*dmats1[5][6] + coeff0_6*dmats1[6][6] + coeff0_7*dmats1[7][6] + coeff0_8*dmats1[8][6] + coeff0_9*dmats1[9][6];
          new_coeff0_7 = coeff0_0*dmats1[0][7] + coeff0_1*dmats1[1][7] + coeff0_2*dmats1[2][7] + coeff0_3*dmats1[3][7] + coeff0_4*dmats1[4][7] + coeff0_5*dmats1[5][7] + coeff0_6*dmats1[6][7] + coeff0_7*dmats1[7][7] + coeff0_8*dmats1[8][7] + coeff0_9*dmats1[9][7];
          new_coeff0_8 = coeff0_0*dmats1[0][8] + coeff0_1*dmats1[1][8] + coeff0_2*dmats1[2][8] + coeff0_3*dmats1[3][8] + coeff0_4*dmats1[4][8] + coeff0_5*dmats1[5][8] + coeff0_6*dmats1[6][8] + coeff0_7*dmats1[7][8] + coeff0_8*dmats1[8][8] + coeff0_9*dmats1[9][8];
          new_coeff0_9 = coeff0_0*dmats1[0][9] + coeff0_1*dmats1[1][9] + coeff0_2*dmats1[2][9] + coeff0_3*dmats1[3][9] + coeff0_4*dmats1[4][9] + coeff0_5*dmats1[5][9] + coeff0_6*dmats1[6][9] + coeff0_7*dmats1[7][9] + coeff0_8*dmats1[8][9] + coeff0_9*dmats1[9][9];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5 + new_coeff0_6*basisvalue6 + new_coeff0_7*basisvalue7 + new_coeff0_8*basisvalue8 + new_coeff0_9*basisvalue9;
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
    double coordinates[2];
    
    // Nodal coordinates on reference cell
    static double X[10][2] = {{0, 0}, {1, 0}, {0, 1}, {0.666666666666667, 0.333333333333333}, {0.333333333333333, 0.666666666666667}, {0, 0.666666666666667}, {0, 0.333333333333333}, {0.333333333333333, 0}, {0.666666666666667, 0}, {0.333333333333333, 0.333333333333333}};
    
    // Components for each dof
    static unsigned int components[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Evaluate basis functions for affine mapping
    const double w0 = 1.0 - X[i][0] - X[i][1];
    const double w1 = X[i][0];
    const double w2 = X[i][1];
    
    // Compute affine mapping x = F(X)
    coordinates[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
    coordinates[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];
    
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
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new ffc_02_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_02_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  ffc_02_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~ffc_02_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Lagrange finite element of degree 3 on a triangle";
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
      return true;
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
    __global_dimension = m.num_entities[0] + 2*m.num_entities[1] + m.num_entities[2];
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
    return 10;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 4;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    unsigned int offset = m.num_entities[0];
    dofs[3] = offset + 2*c.entity_indices[1][0];
    dofs[4] = offset + 2*c.entity_indices[1][0] + 1;
    dofs[5] = offset + 2*c.entity_indices[1][1];
    dofs[6] = offset + 2*c.entity_indices[1][1] + 1;
    dofs[7] = offset + 2*c.entity_indices[1][2];
    dofs[8] = offset + 2*c.entity_indices[1][2] + 1;
    offset = offset + 2*m.num_entities[1];
    dofs[9] = offset + c.entity_indices[2][0];
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   unsigned int facet) const
  {
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 5;
      dofs[3] = 6;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 7;
      dofs[3] = 8;
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::cell& c) const
  {
    const double * const * x = c.coordinates;
    coordinates[0][0] = x[0][0];
    coordinates[0][1] = x[0][1];
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[3][0] = 0.666666666666667*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[3][1] = 0.666666666666667*x[1][1] + 0.333333333333333*x[2][1];
    coordinates[4][0] = 0.333333333333333*x[1][0] + 0.666666666666667*x[2][0];
    coordinates[4][1] = 0.333333333333333*x[1][1] + 0.666666666666667*x[2][1];
    coordinates[5][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[2][0];
    coordinates[5][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[2][1];
    coordinates[6][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[2][0];
    coordinates[6][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[2][1];
    coordinates[7][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[1][0];
    coordinates[7][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[1][1];
    coordinates[8][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[1][0];
    coordinates[8][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[1][1];
    coordinates[9][0] = 0.333333333333333*x[0][0] + 0.333333333333333*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[9][1] = 0.333333333333333*x[0][1] + 0.333333333333333*x[1][1] + 0.333333333333333*x[2][1];
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new ffc_02_dof_map_0();
  }

};

#endif
