// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.7.0.
//
// Warning: This code was generated with the option '-l dolfin'
// and contains DOLFIN-specific wrappers that depend on DOLFIN.

#ifndef __P3_H
#define __P3_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>
    
/// This class defines the interface for a finite element.

class p3_0_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  p3_0_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~p3_0_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "FiniteElement('Lagrange', 'triangle', 3)";
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
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
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
    static const double coefficients0[10][10] = \
    {{0.047140452079103, -0.0288675134594813, -0.0166666666666666, 0.0782460796435951, 0.0606091526731326, 0.0349927106111883, -0.0601337794302955, -0.0508223195384204, -0.0393667994375868, -0.0227284322524248},
    {0.0471404520791031, 0.0288675134594813, -0.0166666666666666, 0.0782460796435952, -0.0606091526731326, 0.0349927106111883, 0.0601337794302955, -0.0508223195384204, 0.0393667994375868, -0.0227284322524248},
    {0.0471404520791032, 0, 0.0333333333333334, 0, 0, 0.104978131833565, 0, 0, 0, 0.0909137290096989},
    {0.106066017177982, 0.259807621135332, -0.15, 0.117369119465393, 0.0606091526731326, -0.0787335988751736, 0, 0.101644639076841, -0.131222664791956, 0.090913729009699},
    {0.106066017177982, 0, 0.3, 0, 0.151522881682832, 0.0262445329583912, 0, 0, 0.131222664791956, -0.136370593514548},
    {0.106066017177982, -0.259807621135332, -0.15, 0.117369119465393, -0.0606091526731326, -0.0787335988751736, 0, 0.101644639076841, 0.131222664791956, 0.090913729009699},
    {0.106066017177982, 0, 0.3, 0, -0.151522881682832, 0.0262445329583912, 0, 0, -0.131222664791956, -0.136370593514548},
    {0.106066017177982, -0.259807621135332, -0.15, -0.0782460796435951, 0.090913729009699, 0.0962299541807677, 0.180401338290886, 0.0508223195384204, -0.0131222664791956, -0.0227284322524247},
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
    const double J_10 = element_coordinates[1][1] - element_coordinates[0][1];
    const double J_11 = element_coordinates[2][1] - element_coordinates[0][1];
    
    // Compute determinant of Jacobian
    const double detJ = J_00*J_11 - J_01*J_10;
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (UFC) element
    double x = (element_coordinates[0][1]*element_coordinates[2][0] -\
                element_coordinates[0][0]*element_coordinates[2][1] +\
                J_11*coordinates[0] - J_01*coordinates[1]) / detJ;
    double y = (element_coordinates[1][1]*element_coordinates[0][0] -\
                element_coordinates[1][0]*element_coordinates[0][1] -\
                J_10*coordinates[0] + J_00*coordinates[1]) / detJ;
    
    // Map coordinates to the reference square
    if (std::abs(y - 1.0) < 1e-14)
      x = -1.0;
    else
      x = 2.0 *x/(1.0 - y) - 1.0;
    y = 2.0*y - 1.0;
    
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
    
    // Compute inverse of Jacobian
    const double Jinv[2][2] =  {{J_11 / detJ, -J_01 / detJ}, {-J_10 / detJ, J_00 / detJ}};
    
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
    static const double coefficients0[10][10] = \
    {{0.047140452079103, -0.0288675134594813, -0.0166666666666666, 0.0782460796435951, 0.0606091526731326, 0.0349927106111883, -0.0601337794302955, -0.0508223195384204, -0.0393667994375868, -0.0227284322524248},
    {0.0471404520791031, 0.0288675134594813, -0.0166666666666666, 0.0782460796435952, -0.0606091526731326, 0.0349927106111883, 0.0601337794302955, -0.0508223195384204, 0.0393667994375868, -0.0227284322524248},
    {0.0471404520791032, 0, 0.0333333333333334, 0, 0, 0.104978131833565, 0, 0, 0, 0.0909137290096989},
    {0.106066017177982, 0.259807621135332, -0.15, 0.117369119465393, 0.0606091526731326, -0.0787335988751736, 0, 0.101644639076841, -0.131222664791956, 0.090913729009699},
    {0.106066017177982, 0, 0.3, 0, 0.151522881682832, 0.0262445329583912, 0, 0, 0.131222664791956, -0.136370593514548},
    {0.106066017177982, -0.259807621135332, -0.15, 0.117369119465393, -0.0606091526731326, -0.0787335988751736, 0, 0.101644639076841, 0.131222664791956, 0.090913729009699},
    {0.106066017177982, 0, 0.3, 0, -0.151522881682832, 0.0262445329583912, 0, 0, -0.131222664791956, -0.136370593514548},
    {0.106066017177982, -0.259807621135332, -0.15, -0.0782460796435951, 0.090913729009699, 0.0962299541807677, 0.180401338290886, 0.0508223195384204, -0.0131222664791956, -0.0227284322524247},
    {0.106066017177982, 0.259807621135332, -0.15, -0.0782460796435952, -0.090913729009699, 0.0962299541807677, -0.180401338290886, 0.0508223195384204, 0.0131222664791956, -0.0227284322524247},
    {0.636396103067893, 0, 0, -0.234738238930785, 0, -0.262445329583912, 0, -0.203289278153682, 0, 0.090913729009699}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {4.89897948556636, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 9.48683298050514, 0, 0, 0, 0, 0, 0, 0, 0},
    {4, 0, 7.07106781186547, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {5.29150262212918, 0, -2.99332590941915, 13.6626010212795, 0, 0.61101009266078, 0, 0, 0, 0},
    {0, 4.38178046004133, 0, 0, 12.5219806739988, 0, 0, 0, 0, 0},
    {3.46410161513775, 0, 7.83836717690617, 0, 0, 8.4, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    static const double dmats1[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.44948974278318, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {4.24264068711929, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.58198889747161, 4.74341649025257, -0.912870929175274, 0, 0, 0, 0, 0, 0, 0},
    {2, 6.12372435695795, 3.53553390593274, 0, 0, 0, 0, 0, 0, 0},
    {-2.30940107675851, 0, 8.16496580927726, 0, 0, 0, 0, 0, 0, 0},
    {2.64575131106459, 5.18459255872629, -1.49666295470957, 6.83130051063973, -1.05830052442584, 0.305505046330391, 0, 0, 0, 0},
    {2.23606797749979, 2.19089023002067, 2.5298221281347, 8.08290376865476, 6.26099033699941, -1.80739222823013, 0, 0, 0, 0},
    {1.73205080756887, -5.09116882454314, 3.91918358845308, 0, 9.69948452238571, 4.2, 0, 0, 0, 0},
    {5.00000000000001, 0, -2.82842712474619, 0, 0, 12.1243556529821, 0, 0, 0, 0}};
    
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
    static const double X[10][1][2] = {{{0, 0}}, {{1, 0}}, {{0, 1}}, {{0.666666666666667, 0.333333333333333}}, {{0.333333333333333, 0.666666666666667}}, {{0, 0.333333333333333}}, {{0, 0.666666666666667}}, {{0.333333333333333, 0}}, {{0.666666666666667, 0}}, {{0.333333333333333, 0.333333333333333}}};
    static const double W[10][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    static const double D[10][1][1] = {{{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}};
    
    const double * const * x = c.coordinates;
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
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new p3_0_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class p3_0_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  p3_0_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~p3_0_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for FiniteElement('Lagrange', 'triangle', 3)";
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

  /// Return the dimension of the local finite element function space for a cell
  virtual unsigned int local_dimension(const ufc::cell& c) const
  {
    return 10;
  }

  /// Return the maximum dimension of the local finite element function space
  virtual unsigned int max_local_dimension() const
  {
    return 10;
  }

  // Return the geometric dimension of the coordinates this dof map provides
  virtual unsigned int geometric_dimension() const
  {
    return 2;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 4;
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
    coordinates[1][0] = x[1][0];
    coordinates[1][1] = x[1][1];
    coordinates[2][0] = x[2][0];
    coordinates[2][1] = x[2][1];
    coordinates[3][0] = 0.666666666666667*x[1][0] + 0.333333333333333*x[2][0];
    coordinates[3][1] = 0.666666666666667*x[1][1] + 0.333333333333333*x[2][1];
    coordinates[4][0] = 0.333333333333333*x[1][0] + 0.666666666666667*x[2][0];
    coordinates[4][1] = 0.333333333333333*x[1][1] + 0.666666666666667*x[2][1];
    coordinates[5][0] = 0.666666666666667*x[0][0] + 0.333333333333333*x[2][0];
    coordinates[5][1] = 0.666666666666667*x[0][1] + 0.333333333333333*x[2][1];
    coordinates[6][0] = 0.333333333333333*x[0][0] + 0.666666666666667*x[2][0];
    coordinates[6][1] = 0.333333333333333*x[0][1] + 0.666666666666667*x[2][1];
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
    return new p3_0_dof_map_0();
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

class p3_form_0: public ufc::form
{
public:

  /// Constructor
  p3_form_0() : ufc::form()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~p3_form_0()
  {
    // Do nothing
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
    return "None";
  }

  /// Return the rank of the global tensor (r)
  virtual unsigned int rank() const
  {
    return -1;
  }

  /// Return the number of coefficients (n)
  virtual unsigned int num_coefficients() const
  {
    return 0;
  }

  /// Return the number of cell integrals
  virtual unsigned int num_cell_integrals() const
  {
    return 0;
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
    return 0;
  }

  /// Create a new dof map for argument function i
  virtual ufc::dof_map* create_dof_map(unsigned int i) const
  {
    return 0;
  }

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(unsigned int i) const
  {
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

namespace P3
{

class Form_0_FunctionSpace_0: public dolfin::FunctionSpace
{
public:

  Form_0_FunctionSpace_0(const dolfin::Mesh& mesh):
      dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new p3_0_finite_element_0()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new p3_0_dof_map_0()), dolfin::reference_to_no_delete_pointer(mesh))))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_0(dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new p3_0_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new p3_0_dof_map_0()), dolfin::reference_to_no_delete_pointer(mesh))))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_0(boost::shared_ptr<dolfin::Mesh> mesh):
      dolfin::FunctionSpace(mesh,
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new p3_0_finite_element_0()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new p3_0_dof_map_0()), mesh)))
  {
      // Do nothing
  }

  Form_0_FunctionSpace_0(boost::shared_ptr<const dolfin::Mesh> mesh):
      dolfin::FunctionSpace(mesh,
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new p3_0_finite_element_0()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new p3_0_dof_map_0()), mesh)))
  {
      // Do nothing
  }


  ~Form_0_FunctionSpace_0()
  {
  }

};

class Form_0: public dolfin::Form
{
public:

  // Constructor
  Form_0(const dolfin::FunctionSpace& V0):
    dolfin::Form(1, 0)
  {
    _function_spaces[0] = reference_to_no_delete_pointer(V0);

    _ufc_form = boost::shared_ptr<const ufc::form>(new p3_form_0());
  }

  // Constructor
  Form_0(boost::shared_ptr<const dolfin::FunctionSpace> V0):
    dolfin::Form(1, 0)
  {
    _function_spaces[0] = V0;

    _ufc_form = boost::shared_ptr<const ufc::form>(new p3_form_0());
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

  // Coefficients
};

// Class typedefs
typedef Form_0 LinearForm;
typedef Form_0::TestSpace FunctionSpace;

} // namespace P3

#endif
