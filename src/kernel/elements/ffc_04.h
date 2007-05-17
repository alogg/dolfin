// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.3.5.

#ifndef __FFC_04_H
#define __FFC_04_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_04_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_04_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_04_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Lagrange finite element of degree 2 on a tetrahedron";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::tetrahedron;
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
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    const double scalings_z_2 = scalings_z_1*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_00_2 = 0.3125*psitilde_cs_00_1 + 1.875*z*psitilde_cs_00_1 - 0.5625*psitilde_cs_00_0;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_01_1 = 3*z + 2;
    const double psitilde_cs_02_0 = 1;
    const double psitilde_cs_10_0 = 1;
    const double psitilde_cs_10_1 = 3*z + 2;
    const double psitilde_cs_11_0 = 1;
    const double psitilde_cs_20_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    const double basisvalue4 = 5.1234753829798*psitilde_a_2*scalings_y_2*psitilde_bs_2_0*scalings_z_2*psitilde_cs_20_0;
    const double basisvalue5 = 3.96862696659689*psitilde_a_1*scalings_y_1*psitilde_bs_1_1*scalings_z_2*psitilde_cs_11_0;
    const double basisvalue6 = 2.29128784747792*psitilde_a_0*scalings_y_0*psitilde_bs_0_2*scalings_z_2*psitilde_cs_02_0;
    const double basisvalue7 = 3.24037034920393*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_1;
    const double basisvalue8 = 1.87082869338697*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_1;
    const double basisvalue9 = 1.3228756555323*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_2;
    
    // Table(s) of coefficients
    const static double coefficients0[10][10] = \
    {{-0.0577350269189626, -0.0608580619450185, -0.0351364184463153, -0.0248451997499977, 0.0650600048632355, 0.050395263067897, 0.0290957186981323, 0.0411475599898912, 0.0237565548366599, 0.0167984210226323},
    {-0.0577350269189625, 0.0608580619450185, -0.0351364184463153, -0.0248451997499977, 0.0650600048632355, -0.050395263067897, 0.0290957186981323, -0.0411475599898912, 0.0237565548366599, 0.0167984210226323},
    {-0.0577350269189626, 0, 0.0702728368926307, -0.0248451997499977, 0, 0, 0.0872871560943969, 0, -0.0475131096733199, 0.0167984210226323},
    {-0.0577350269189626, 0, 0, 0.074535599249993, 0, 0, 0, 0, 0, 0.100790526135794},
    {0.23094010767585, 0.121716123890037, 0.0702728368926306, -0.0993807989999906, 0, 0.100790526135794, -0.0872871560943969, -0.0205737799949456, -0.01187827741833, 0.0167984210226323},
    {0.23094010767585, -0.121716123890037, 0.0702728368926306, -0.0993807989999907, 0, -0.100790526135794, -0.0872871560943969, 0.0205737799949456, -0.01187827741833, 0.0167984210226323},
    {0.23094010767585, 0, -0.140545673785261, -0.0993807989999906, -0.130120009726471, 0, 0.0290957186981323, 0, 0.02375655483666, 0.0167984210226323},
    {0.23094010767585, -0.121716123890037, -0.0702728368926306, 0.0993807989999906, 0, 0, 0, -0.102868899974728, -0.0593913870916499, -0.0671936840905293},
    {0.23094010767585, 0.121716123890037, -0.0702728368926307, 0.0993807989999907, 0, 0, 0, 0.102868899974728, -0.0593913870916499, -0.0671936840905293},
    {0.23094010767585, 0, 0.140545673785261, 0.0993807989999906, 0, 0, 0, 0, 0.1187827741833, -0.0671936840905293}};
    
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
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    const double scalings_y_2 = scalings_y_1*(0.5 - 0.5*y);
    const double scalings_z_0 = 1;
    const double scalings_z_1 = scalings_z_0*(0.5 - 0.5*z);
    const double scalings_z_2 = scalings_z_1*(0.5 - 0.5*z);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    const double psitilde_a_2 = 1.5*x*psitilde_a_1 - 0.5*psitilde_a_0;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_0_2 = 0.111111111111111*psitilde_bs_0_1 + 1.66666666666667*y*psitilde_bs_0_1 - 0.555555555555556*psitilde_bs_0_0;
    const double psitilde_bs_1_0 = 1;
    const double psitilde_bs_1_1 = 2.5*y + 1.5;
    const double psitilde_bs_2_0 = 1;
    
    // Compute psitilde_cs
    const double psitilde_cs_00_0 = 1;
    const double psitilde_cs_00_1 = 2*z + 1;
    const double psitilde_cs_00_2 = 0.3125*psitilde_cs_00_1 + 1.875*z*psitilde_cs_00_1 - 0.5625*psitilde_cs_00_0;
    const double psitilde_cs_01_0 = 1;
    const double psitilde_cs_01_1 = 3*z + 2;
    const double psitilde_cs_02_0 = 1;
    const double psitilde_cs_10_0 = 1;
    const double psitilde_cs_10_1 = 3*z + 2;
    const double psitilde_cs_11_0 = 1;
    const double psitilde_cs_20_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.866025403784439*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_0;
    const double basisvalue1 = 2.73861278752583*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_0;
    const double basisvalue2 = 1.58113883008419*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_0;
    const double basisvalue3 = 1.11803398874989*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_1;
    const double basisvalue4 = 5.1234753829798*psitilde_a_2*scalings_y_2*psitilde_bs_2_0*scalings_z_2*psitilde_cs_20_0;
    const double basisvalue5 = 3.96862696659689*psitilde_a_1*scalings_y_1*psitilde_bs_1_1*scalings_z_2*psitilde_cs_11_0;
    const double basisvalue6 = 2.29128784747792*psitilde_a_0*scalings_y_0*psitilde_bs_0_2*scalings_z_2*psitilde_cs_02_0;
    const double basisvalue7 = 3.24037034920393*psitilde_a_1*scalings_y_1*psitilde_bs_1_0*scalings_z_1*psitilde_cs_10_1;
    const double basisvalue8 = 1.87082869338697*psitilde_a_0*scalings_y_0*psitilde_bs_0_1*scalings_z_1*psitilde_cs_01_1;
    const double basisvalue9 = 1.3228756555323*psitilde_a_0*scalings_y_0*psitilde_bs_0_0*scalings_z_0*psitilde_cs_00_2;
    
    // Table(s) of coefficients
    const static double coefficients0[10][10] = \
    {{-0.0577350269189626, -0.0608580619450185, -0.0351364184463153, -0.0248451997499977, 0.0650600048632355, 0.050395263067897, 0.0290957186981323, 0.0411475599898912, 0.0237565548366599, 0.0167984210226323},
    {-0.0577350269189625, 0.0608580619450185, -0.0351364184463153, -0.0248451997499977, 0.0650600048632355, -0.050395263067897, 0.0290957186981323, -0.0411475599898912, 0.0237565548366599, 0.0167984210226323},
    {-0.0577350269189626, 0, 0.0702728368926307, -0.0248451997499977, 0, 0, 0.0872871560943969, 0, -0.0475131096733199, 0.0167984210226323},
    {-0.0577350269189626, 0, 0, 0.074535599249993, 0, 0, 0, 0, 0, 0.100790526135794},
    {0.23094010767585, 0.121716123890037, 0.0702728368926306, -0.0993807989999906, 0, 0.100790526135794, -0.0872871560943969, -0.0205737799949456, -0.01187827741833, 0.0167984210226323},
    {0.23094010767585, -0.121716123890037, 0.0702728368926306, -0.0993807989999907, 0, -0.100790526135794, -0.0872871560943969, 0.0205737799949456, -0.01187827741833, 0.0167984210226323},
    {0.23094010767585, 0, -0.140545673785261, -0.0993807989999906, -0.130120009726471, 0, 0.0290957186981323, 0, 0.02375655483666, 0.0167984210226323},
    {0.23094010767585, -0.121716123890037, -0.0702728368926306, 0.0993807989999906, 0, 0, 0, -0.102868899974728, -0.0593913870916499, -0.0671936840905293},
    {0.23094010767585, 0.121716123890037, -0.0702728368926307, 0.0993807989999907, 0, 0, 0, 0.102868899974728, -0.0593913870916499, -0.0671936840905293},
    {0.23094010767585, 0, 0.140545673785261, 0.0993807989999906, 0, 0, 0, 0, 0.1187827741833, -0.0671936840905293}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {3.16227766016838, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 5.61248608016091, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.29128784747792, 0, 4.18330013267038, -0.591607978309962, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.87082869338697, 0, 0, 4.34741302385683, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    const static double dmats1[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.58113883008419, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.73861278752583, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.4790199457749, 2.80624304008046, -0.540061724867322, -0.381881307912987, 0, 0, 0, 0, 0, 0},
    {1.14564392373896, 3.62284418654736, 2.09165006633519, -0.295803989154981, 0, 0, 0, 0, 0, 0},
    {-1.3228756555323, 0, 4.83045891539648, 0.341565025531987, 0, 0, 0, 0, 0, 0},
    {0.935414346693486, 0, 0, 2.17370651192842, 0, 0, 0, 0, 0, 0},
    {1.62018517460197, 0, 0, 3.76497011940334, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    
    const static double dmats2[10][10] = \
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.58113883008419, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0.912870929175277, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {2.58198889747161, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1.4790199457749, 2.80624304008046, -0.540061724867322, -0.381881307912987, 0, 0, 0, 0, 0, 0},
    {1.14564392373896, 0.724568837309472, 2.09165006633519, -0.295803989154981, 0, 0, 0, 0, 0, 0},
    {0.661437827766147, 0, 1.93218356615859, -0.170782512765993, 0, 0, 0, 0, 0, 0},
    {0.935414346693486, 3.54964786985977, 0, 2.17370651192842, 0, 0, 0, 0, 0, 0},
    {0.540061724867322, 0, 3.54964786985977, 1.25499003980111, 0, 0, 0, 0, 0, 0},
    {-1.90940653956493, 0, 0, 4.43705983732471, 0, 0, 0, 0, 0, 0}};
    
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
        if(combinations[deriv_num][j] == 2)
        {
          new_coeff0_0 = coeff0_0*dmats2[0][0] + coeff0_1*dmats2[1][0] + coeff0_2*dmats2[2][0] + coeff0_3*dmats2[3][0] + coeff0_4*dmats2[4][0] + coeff0_5*dmats2[5][0] + coeff0_6*dmats2[6][0] + coeff0_7*dmats2[7][0] + coeff0_8*dmats2[8][0] + coeff0_9*dmats2[9][0];
          new_coeff0_1 = coeff0_0*dmats2[0][1] + coeff0_1*dmats2[1][1] + coeff0_2*dmats2[2][1] + coeff0_3*dmats2[3][1] + coeff0_4*dmats2[4][1] + coeff0_5*dmats2[5][1] + coeff0_6*dmats2[6][1] + coeff0_7*dmats2[7][1] + coeff0_8*dmats2[8][1] + coeff0_9*dmats2[9][1];
          new_coeff0_2 = coeff0_0*dmats2[0][2] + coeff0_1*dmats2[1][2] + coeff0_2*dmats2[2][2] + coeff0_3*dmats2[3][2] + coeff0_4*dmats2[4][2] + coeff0_5*dmats2[5][2] + coeff0_6*dmats2[6][2] + coeff0_7*dmats2[7][2] + coeff0_8*dmats2[8][2] + coeff0_9*dmats2[9][2];
          new_coeff0_3 = coeff0_0*dmats2[0][3] + coeff0_1*dmats2[1][3] + coeff0_2*dmats2[2][3] + coeff0_3*dmats2[3][3] + coeff0_4*dmats2[4][3] + coeff0_5*dmats2[5][3] + coeff0_6*dmats2[6][3] + coeff0_7*dmats2[7][3] + coeff0_8*dmats2[8][3] + coeff0_9*dmats2[9][3];
          new_coeff0_4 = coeff0_0*dmats2[0][4] + coeff0_1*dmats2[1][4] + coeff0_2*dmats2[2][4] + coeff0_3*dmats2[3][4] + coeff0_4*dmats2[4][4] + coeff0_5*dmats2[5][4] + coeff0_6*dmats2[6][4] + coeff0_7*dmats2[7][4] + coeff0_8*dmats2[8][4] + coeff0_9*dmats2[9][4];
          new_coeff0_5 = coeff0_0*dmats2[0][5] + coeff0_1*dmats2[1][5] + coeff0_2*dmats2[2][5] + coeff0_3*dmats2[3][5] + coeff0_4*dmats2[4][5] + coeff0_5*dmats2[5][5] + coeff0_6*dmats2[6][5] + coeff0_7*dmats2[7][5] + coeff0_8*dmats2[8][5] + coeff0_9*dmats2[9][5];
          new_coeff0_6 = coeff0_0*dmats2[0][6] + coeff0_1*dmats2[1][6] + coeff0_2*dmats2[2][6] + coeff0_3*dmats2[3][6] + coeff0_4*dmats2[4][6] + coeff0_5*dmats2[5][6] + coeff0_6*dmats2[6][6] + coeff0_7*dmats2[7][6] + coeff0_8*dmats2[8][6] + coeff0_9*dmats2[9][6];
          new_coeff0_7 = coeff0_0*dmats2[0][7] + coeff0_1*dmats2[1][7] + coeff0_2*dmats2[2][7] + coeff0_3*dmats2[3][7] + coeff0_4*dmats2[4][7] + coeff0_5*dmats2[5][7] + coeff0_6*dmats2[6][7] + coeff0_7*dmats2[7][7] + coeff0_8*dmats2[8][7] + coeff0_9*dmats2[9][7];
          new_coeff0_8 = coeff0_0*dmats2[0][8] + coeff0_1*dmats2[1][8] + coeff0_2*dmats2[2][8] + coeff0_3*dmats2[3][8] + coeff0_4*dmats2[4][8] + coeff0_5*dmats2[5][8] + coeff0_6*dmats2[6][8] + coeff0_7*dmats2[7][8] + coeff0_8*dmats2[8][8] + coeff0_9*dmats2[9][8];
          new_coeff0_9 = coeff0_0*dmats2[0][9] + coeff0_1*dmats2[1][9] + coeff0_2*dmats2[2][9] + coeff0_3*dmats2[3][9] + coeff0_4*dmats2[4][9] + coeff0_5*dmats2[5][9] + coeff0_6*dmats2[6][9] + coeff0_7*dmats2[7][9] + coeff0_8*dmats2[8][9] + coeff0_9*dmats2[9][9];
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
    double coordinates[3];
    
    // Nodal coordinates on reference cell
    static double X[10][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.5, 0.5, 0}, {0, 0.5, 0}, {0.5, 0, 0}, {0, 0, 0.5}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    
    // Components for each dof
    static unsigned int components[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
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
    return new ffc_04_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_04_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  ffc_04_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~ffc_04_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Lagrange finite element of degree 2 on a tetrahedron";
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
    __global_dimension = m.num_entities[0] + m.num_entities[1];
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
    return 6;
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
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[1][0];
    dofs[5] = offset + c.entity_indices[1][1];
    dofs[6] = offset + c.entity_indices[1][2];
    dofs[7] = offset + c.entity_indices[1][3];
    dofs[8] = offset + c.entity_indices[1][4];
    dofs[9] = offset + c.entity_indices[1][5];
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
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 7;
      dofs[5] = 8;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 7;
      dofs[5] = 9;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 6;
      dofs[4] = 8;
      dofs[5] = 9;
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::cell& c) const
  {
    // This function is implemented assuming affine mapping!!
    // Get cell vertices
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
    coordinates[4][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[4][2] = 0.5*x[1][2] + 0.5*x[2][2];
    coordinates[5][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[5][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][2] = 0.5*x[0][2] + 0.5*x[2][2];
    coordinates[6][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[6][1] = 0.5*x[0][1] + 0.5*x[1][1];
    coordinates[6][2] = 0.5*x[0][2] + 0.5*x[1][2];
    coordinates[7][0] = 0.5*x[0][0] + 0.5*x[3][0];
    coordinates[7][1] = 0.5*x[0][1] + 0.5*x[3][1];
    coordinates[7][2] = 0.5*x[0][2] + 0.5*x[3][2];
    coordinates[8][0] = 0.5*x[1][0] + 0.5*x[3][0];
    coordinates[8][1] = 0.5*x[1][1] + 0.5*x[3][1];
    coordinates[8][2] = 0.5*x[1][2] + 0.5*x[3][2];
    coordinates[9][0] = 0.5*x[2][0] + 0.5*x[3][0];
    coordinates[9][1] = 0.5*x[2][1] + 0.5*x[3][1];
    coordinates[9][2] = 0.5*x[2][2] + 0.5*x[3][2];
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new ffc_04_dof_map_0();
  }

};

#endif
