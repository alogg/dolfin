#include "CahnHilliard3D.h"

/// Constructor
cahnhilliard3d_0_finite_element_0_0::cahnhilliard3d_0_finite_element_0_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_0_0::~cahnhilliard3d_0_finite_element_0_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_0_0::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_0_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_0_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_0_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_0_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_0_0::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_0_finite_element_0_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_0_0::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_0_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_0_0::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_0_finite_element_0_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_0_0::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_0_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_0_0::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_0_0();
}


/// Constructor
cahnhilliard3d_0_finite_element_0_1::cahnhilliard3d_0_finite_element_0_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_0_1::~cahnhilliard3d_0_finite_element_0_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_0_1::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_0_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_0_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_0_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_0_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_0_1::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_0_finite_element_0_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_0_1::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_0_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_0_1::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_0_finite_element_0_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_0_1::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_0_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_0_1::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_0_1();
}


/// Constructor
cahnhilliard3d_0_finite_element_0::cahnhilliard3d_0_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_0::~cahnhilliard3d_0_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_0::signature() const
{
    return "MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_0::space_dimension() const
{
    return 8;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_0::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_0::evaluate_basis(unsigned int i,
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
    values[0] = 0;
    values[1] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_0::evaluate_basis_derivatives(unsigned int i,
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
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
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[8][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[8][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[8][1][2] = {{{1, 0}}, {{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}, {{0, 1}}, {{0, 1}}, {{0, 1}}};
    
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
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void cahnhilliard3d_0_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[2] = dof_values[1];
    vertex_values[4] = dof_values[2];
    vertex_values[6] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[3] = dof_values[5];
    vertex_values[5] = dof_values[6];
    vertex_values[7] = dof_values[7];
}

/// Return the number of sub elements (for a mixed element)
unsigned int cahnhilliard3d_0_finite_element_0::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_0::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_finite_element_0_0();
      break;
    case 1:
      return new cahnhilliard3d_0_finite_element_0_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_0_finite_element_1_0::cahnhilliard3d_0_finite_element_1_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_1_0::~cahnhilliard3d_0_finite_element_1_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_1_0::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_1_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_1_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_1_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_1_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_1_0::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_0_finite_element_1_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_1_0::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_1_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_1_0::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_0_finite_element_1_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_1_0::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_1_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_1_0::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_1_0();
}


/// Constructor
cahnhilliard3d_0_finite_element_1_1::cahnhilliard3d_0_finite_element_1_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_1_1::~cahnhilliard3d_0_finite_element_1_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_1_1::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_1_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_1_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_1_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_1_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_1_1::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_0_finite_element_1_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_1_1::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_1_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_1_1::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_0_finite_element_1_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_1_1::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_1_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_1_1::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_1_1();
}


/// Constructor
cahnhilliard3d_0_finite_element_1::cahnhilliard3d_0_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_1::~cahnhilliard3d_0_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_1::signature() const
{
    return "MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_1::space_dimension() const
{
    return 8;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_1::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_1::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_1::evaluate_basis(unsigned int i,
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
    values[0] = 0;
    values[1] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_1::evaluate_basis_derivatives(unsigned int i,
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
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
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[8][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[8][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[8][1][2] = {{{1, 0}}, {{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}, {{0, 1}}, {{0, 1}}, {{0, 1}}};
    
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
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void cahnhilliard3d_0_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[2] = dof_values[1];
    vertex_values[4] = dof_values[2];
    vertex_values[6] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[3] = dof_values[5];
    vertex_values[5] = dof_values[6];
    vertex_values[7] = dof_values[7];
}

/// Return the number of sub elements (for a mixed element)
unsigned int cahnhilliard3d_0_finite_element_1::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_1::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_finite_element_1_0();
      break;
    case 1:
      return new cahnhilliard3d_0_finite_element_1_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_0_finite_element_2_0::cahnhilliard3d_0_finite_element_2_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_2_0::~cahnhilliard3d_0_finite_element_2_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_2_0::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_2_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_2_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_2_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_2_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_2_0::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_0_finite_element_2_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_2_0::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_2_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_2_0::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_0_finite_element_2_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_2_0::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_2_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_2_0::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_2_0();
}


/// Constructor
cahnhilliard3d_0_finite_element_2_1::cahnhilliard3d_0_finite_element_2_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_2_1::~cahnhilliard3d_0_finite_element_2_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_2_1::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_2_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_2_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_2_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_2_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_2_1::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_0_finite_element_2_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_2_1::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_2_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_2_1::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_0_finite_element_2_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_2_1::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_2_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_2_1::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_2_1();
}


/// Constructor
cahnhilliard3d_0_finite_element_2::cahnhilliard3d_0_finite_element_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_2::~cahnhilliard3d_0_finite_element_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_2::signature() const
{
    return "MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_2::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_2::space_dimension() const
{
    return 8;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_2::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_2::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_2::evaluate_basis(unsigned int i,
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
    values[0] = 0;
    values[1] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_2::evaluate_basis_derivatives(unsigned int i,
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
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
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_2::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[8][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[8][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[8][1][2] = {{{1, 0}}, {{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}, {{0, 1}}, {{0, 1}}, {{0, 1}}};
    
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
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void cahnhilliard3d_0_finite_element_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[2] = dof_values[1];
    vertex_values[4] = dof_values[2];
    vertex_values[6] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[3] = dof_values[5];
    vertex_values[5] = dof_values[6];
    vertex_values[7] = dof_values[7];
}

/// Return the number of sub elements (for a mixed element)
unsigned int cahnhilliard3d_0_finite_element_2::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_2::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_finite_element_2_0();
      break;
    case 1:
      return new cahnhilliard3d_0_finite_element_2_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_0_finite_element_3::cahnhilliard3d_0_finite_element_3() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_3::~cahnhilliard3d_0_finite_element_3()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_3::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_3::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_3::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_3::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_3::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_3::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_3::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_3::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_3::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_3::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_0_finite_element_3::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_3::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_3::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_3::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_3();
}


/// Constructor
cahnhilliard3d_0_finite_element_4::cahnhilliard3d_0_finite_element_4() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_4::~cahnhilliard3d_0_finite_element_4()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_4::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_4::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_4::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_4::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_4::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_4::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_4::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_4::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_4::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_4::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_0_finite_element_4::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_4::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_4::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_4::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_4();
}


/// Constructor
cahnhilliard3d_0_finite_element_5::cahnhilliard3d_0_finite_element_5() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_5::~cahnhilliard3d_0_finite_element_5()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_5::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_5::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_5::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_5::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_5::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_5::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_5::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_5::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_5::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_5::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_0_finite_element_5::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_5::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_5::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_5::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_5();
}


/// Constructor
cahnhilliard3d_0_finite_element_6::cahnhilliard3d_0_finite_element_6() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_finite_element_6::~cahnhilliard3d_0_finite_element_6()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_0_finite_element_6::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_0_finite_element_6::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_0_finite_element_6::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_0_finite_element_6::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_0_finite_element_6::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_0_finite_element_6::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_0_finite_element_6::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_0_finite_element_6::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_0_finite_element_6::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_0_finite_element_6::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_0_finite_element_6::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_0_finite_element_6::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_0_finite_element_6::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_0_finite_element_6::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_0_finite_element_6();
}

/// Constructor
cahnhilliard3d_0_dof_map_0_0::cahnhilliard3d_0_dof_map_0_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_0_0::~cahnhilliard3d_0_dof_map_0_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_0_0::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_0_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_0_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_0_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_0_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_0_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_0_0::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_0_0::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_0_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_0_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_0_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_0_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_0_0::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_0_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_0_0::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_0_dof_map_0_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_0_0::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_0_0();
}


/// Constructor
cahnhilliard3d_0_dof_map_0_1::cahnhilliard3d_0_dof_map_0_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_0_1::~cahnhilliard3d_0_dof_map_0_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_0_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_0_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_0_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_0_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_0_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_0_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_0_1::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_0_1::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_0_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_0_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_0_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_0_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_0_1::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_0_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_0_1::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_0_dof_map_0_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_0_1::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_0_1();
}


/// Constructor
cahnhilliard3d_0_dof_map_0::cahnhilliard3d_0_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_0::~cahnhilliard3d_0_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_0::signature() const
{
    return "FFC dof map for MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_0::local_dimension(const ufc::cell& c) const
{
    return 8;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_0::max_local_dimension() const
{
    return 8;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_0::num_facet_dofs() const
{
    return 6;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_0::tabulate_coordinates(double** coordinates,
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
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_0::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_dof_map_0_0();
      break;
    case 1:
      return new cahnhilliard3d_0_dof_map_0_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_0_dof_map_1_0::cahnhilliard3d_0_dof_map_1_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_1_0::~cahnhilliard3d_0_dof_map_1_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_1_0::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_1_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_1_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_1_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_1_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_1_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_1_0::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_1_0::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_1_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_1_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_1_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_1_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_1_0::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_1_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_1_0::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_0_dof_map_1_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_1_0::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_1_0();
}


/// Constructor
cahnhilliard3d_0_dof_map_1_1::cahnhilliard3d_0_dof_map_1_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_1_1::~cahnhilliard3d_0_dof_map_1_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_1_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_1_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_1_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_1_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_1_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_1_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_1_1::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_1_1::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_1_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_1_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_1_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_1_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_1_1::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_1_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_1_1::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_0_dof_map_1_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_1_1::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_1_1();
}


/// Constructor
cahnhilliard3d_0_dof_map_1::cahnhilliard3d_0_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_1::~cahnhilliard3d_0_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_1::signature() const
{
    return "FFC dof map for MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_1::local_dimension(const ufc::cell& c) const
{
    return 8;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_1::max_local_dimension() const
{
    return 8;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_1::num_facet_dofs() const
{
    return 6;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_1::tabulate_coordinates(double** coordinates,
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
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_1::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_dof_map_1_0();
      break;
    case 1:
      return new cahnhilliard3d_0_dof_map_1_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_0_dof_map_2_0::cahnhilliard3d_0_dof_map_2_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_2_0::~cahnhilliard3d_0_dof_map_2_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_2_0::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_2_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_2_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_2_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_2_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_2_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_2_0::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_2_0::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_2_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_2_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_2_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_2_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_2_0::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_2_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_2_0::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_0_dof_map_2_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_2_0::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_2_0();
}


/// Constructor
cahnhilliard3d_0_dof_map_2_1::cahnhilliard3d_0_dof_map_2_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_2_1::~cahnhilliard3d_0_dof_map_2_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_2_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_2_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_2_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_2_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_2_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_2_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_2_1::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_2_1::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_2_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_2_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_2_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_2_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_2_1::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_2_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_2_1::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_0_dof_map_2_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_2_1::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_2_1();
}


/// Constructor
cahnhilliard3d_0_dof_map_2::cahnhilliard3d_0_dof_map_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_2::~cahnhilliard3d_0_dof_map_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_2::signature() const
{
    return "FFC dof map for MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_2::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_2::local_dimension(const ufc::cell& c) const
{
    return 8;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_2::max_local_dimension() const
{
    return 8;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_2::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_2::num_facet_dofs() const
{
    return 6;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_2::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_0_dof_map_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_2::tabulate_coordinates(double** coordinates,
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
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_2::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_2::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_dof_map_2_0();
      break;
    case 1:
      return new cahnhilliard3d_0_dof_map_2_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_0_dof_map_3::cahnhilliard3d_0_dof_map_3() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_3::~cahnhilliard3d_0_dof_map_3()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_3::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_3::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_3::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_3::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_3::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_3::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_3::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_3::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_3::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_3::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_3::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_3::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_3::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_0_dof_map_3::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_3::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_3::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_3::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_3();
}


/// Constructor
cahnhilliard3d_0_dof_map_4::cahnhilliard3d_0_dof_map_4() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_4::~cahnhilliard3d_0_dof_map_4()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_4::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_4::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_4::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_4::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_4::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_4::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_4::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_4::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_4::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_4::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_4::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_4::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_4::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_0_dof_map_4::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_4::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_4::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_4::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_4();
}


/// Constructor
cahnhilliard3d_0_dof_map_5::cahnhilliard3d_0_dof_map_5() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_5::~cahnhilliard3d_0_dof_map_5()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_5::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_5::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_5::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_5::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_5::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_5::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_5::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_5::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_5::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_5::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_5::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_5::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_5::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_0_dof_map_5::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_5::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_5::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_5::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_5();
}


/// Constructor
cahnhilliard3d_0_dof_map_6::cahnhilliard3d_0_dof_map_6() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_0_dof_map_6::~cahnhilliard3d_0_dof_map_6()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_0_dof_map_6::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_0_dof_map_6::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_0_dof_map_6::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_0_dof_map_6::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_0_dof_map_6::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_0_dof_map_6::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_0_dof_map_6::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_0_dof_map_6::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_0_dof_map_6::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_0_dof_map_6::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_0_dof_map_6::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_0_dof_map_6::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_0_dof_map_6::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_0_dof_map_6::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_0_dof_map_6::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_0_dof_map_6::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_0_dof_map_6::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_0_dof_map_6();
}


/// Constructor
cahnhilliard3d_0_cell_integral_0_quadrature::cahnhilliard3d_0_cell_integral_0_quadrature() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_cell_integral_0_quadrature::~cahnhilliard3d_0_cell_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void cahnhilliard3d_0_cell_integral_0_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
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
    const double Jinv_00 = d_00 / detJ;
    const double Jinv_01 = d_10 / detJ;
    const double Jinv_02 = d_20 / detJ;
    const double Jinv_10 = d_01 / detJ;
    const double Jinv_11 = d_11 / detJ;
    const double Jinv_12 = d_21 / detJ;
    const double Jinv_20 = d_02 / detJ;
    const double Jinv_21 = d_12 / detJ;
    const double Jinv_22 = d_22 / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    
    // Array of quadrature weights
    const static double W27[27] = {0.00877047492965105, 0.00816265076654668, 0.0016716811314837, 0.0100061425721761, 0.00931268237947044, 0.00190720341498178, 0.00304787709051818, 0.00283664869563092, 0.000580935315837384, 0.0140327598874417, 0.0130602412264747, 0.00267468981037392, 0.0160098281154818, 0.0149002918071527, 0.00305152546397085, 0.0048766033448291, 0.00453863791300947, 0.000929496505339815, 0.00877047492965105, 0.00816265076654668, 0.0016716811314837, 0.0100061425721761, 0.00931268237947044, 0.00190720341498178, 0.00304787709051818, 0.00283664869563092, 0.000580935315837384};
    
    
    const static double FE0_C1_D001[27][2] = \
    {{-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1}};
    // Array of non-zero columns
    static const unsigned int nzc0[2] = {4, 7};
    // Array of non-zero columns
    static const unsigned int nzc1[2] = {4, 6};
    // Array of non-zero columns
    static const unsigned int nzc2[2] = {4, 5};
    // Array of non-zero columns
    static const unsigned int nzc3[2] = {0, 1};
    // Array of non-zero columns
    static const unsigned int nzc4[2] = {0, 2};
    // Array of non-zero columns
    static const unsigned int nzc5[2] = {0, 3};
    
    const static double FE0_C1[27][4] = \
    {{0.749664528221693, 0.0952198798417149, 0.0821215678634426, 0.0729940240731498},
    {0.528074388273447, 0.0670742417520586, 0.0578476039361428, 0.347003766038352},
    {0.238563056650491, 0.0303014811742758, 0.0261332522867349, 0.705002209888498},
    {0.485731727037113, 0.0616960186091465, 0.379578230280591, 0.0729940240731498},
    {0.342156357895961, 0.0434595556538024, 0.267380320411885, 0.347003766038352},
    {0.154572667042115, 0.0196333029354846, 0.120791820133903, 0.705002209888498},
    {0.174656645238399, 0.0221843026408197, 0.730165028047632, 0.0729940240731498},
    {0.123030632529655, 0.0156269392579017, 0.514338662174092, 0.347003766038352},
    {0.0555803583920822, 0.00705963113955482, 0.232357800579865, 0.705002209888498},
    {0.422442204031704, 0.422442204031704, 0.0821215678634426, 0.0729940240731498},
    {0.297574315012753, 0.297574315012753, 0.0578476039361428, 0.347003766038352},
    {0.134432268912383, 0.134432268912383, 0.0261332522867349, 0.705002209888498},
    {0.27371387282313, 0.27371387282313, 0.379578230280591, 0.0729940240731498},
    {0.192807956774882, 0.192807956774882, 0.267380320411885, 0.347003766038352},
    {0.0871029849887995, 0.0871029849887996, 0.120791820133903, 0.705002209888498},
    {0.0984204739396093, 0.0984204739396094, 0.730165028047632, 0.0729940240731498},
    {0.0693287858937781, 0.0693287858937782, 0.514338662174092, 0.347003766038352},
    {0.0313199947658185, 0.0313199947658185, 0.232357800579865, 0.705002209888498},
    {0.0952198798417148, 0.749664528221693, 0.0821215678634426, 0.0729940240731498},
    {0.0670742417520584, 0.528074388273447, 0.0578476039361428, 0.347003766038352},
    {0.0303014811742757, 0.238563056650491, 0.0261332522867349, 0.705002209888498},
    {0.0616960186091464, 0.485731727037113, 0.379578230280591, 0.0729940240731498},
    {0.0434595556538024, 0.342156357895961, 0.267380320411885, 0.347003766038352},
    {0.0196333029354845, 0.154572667042115, 0.120791820133903, 0.705002209888498},
    {0.0221843026408197, 0.174656645238399, 0.730165028047632, 0.0729940240731498},
    {0.0156269392579016, 0.123030632529655, 0.514338662174092, 0.347003766038352},
    {0.00705963113955483, 0.0555803583920822, 0.232357800579865, 0.705002209888498}};
    // Array of non-zero columns
    static const unsigned int nzc6[4] = {4, 5, 6, 7};
    // Array of non-zero columns
    static const unsigned int nzc7[4] = {0, 1, 2, 3};
    
    // Number of operations to compute geometry constants: 104
    const double G0 =  - det*w[1][0]*(Jinv_10*Jinv_10 + Jinv_11*Jinv_11 + Jinv_12*Jinv_12);
    const double G1 = det*w[3][0]*w[4][0]*(Jinv_00*Jinv_10 + Jinv_01*Jinv_11 + Jinv_02*Jinv_12);
    const double G2 =  - det*w[1][0]*(Jinv_00*Jinv_20 + Jinv_01*Jinv_21 + Jinv_02*Jinv_22);
    const double G3 =  - det*w[1][0]*(Jinv_00*Jinv_10 + Jinv_01*Jinv_11 + Jinv_02*Jinv_12);
    const double G4 =  - 2*det*w[2][0];
    const double G5 = 12*det*w[2][0];
    const double G6 =  - 12*det*w[2][0];
    const double G7 = det*w[3][0]*w[4][0]*(Jinv_00*Jinv_20 + Jinv_01*Jinv_21 + Jinv_02*Jinv_22);
    const double G8 = det*w[3][0]*w[4][0]*(Jinv_10*Jinv_10 + Jinv_11*Jinv_11 + Jinv_12*Jinv_12);
    const double G9 = det*w[3][0]*w[4][0]*(Jinv_10*Jinv_20 + Jinv_11*Jinv_21 + Jinv_12*Jinv_22);
    const double G10 =  - det*w[1][0]*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01 + Jinv_02*Jinv_02);
    const double G11 =  - det*w[1][0]*(Jinv_10*Jinv_20 + Jinv_11*Jinv_21 + Jinv_12*Jinv_22);
    const double G12 = det*w[3][0]*w[4][0]*(Jinv_20*Jinv_20 + Jinv_21*Jinv_21 + Jinv_22*Jinv_22);
    const double G13 =  - det*w[1][0]*(Jinv_20*Jinv_20 + Jinv_21*Jinv_21 + Jinv_22*Jinv_22);
    const double G14 = det*w[3][0]*w[4][0]*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01 + Jinv_02*Jinv_02);
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', True), ('ignore zero tables', True), ('non zero columns', True), ('remove zero terms', True), ('ignore ones', True)
    // Total number of operations to compute element tensor: 10526
    
    // Loop quadrature points for integral
    // Number of operations to compute element tensor for following IP loop = 10422
    for (unsigned int ip = 0; ip < 27; ip++)
    {
      
      // Function declarations
      double F0 = 0;
      
      // Total number of operations to compute function values = 8
      for (unsigned int r = 0; r < 4; r++)
      {
        F0 += FE0_C1[ip][r]*w[0][nzc6[r]];
      }// end loop over 'r'
      
      // Number of operations to compute ip constants: 18
      // Number of operations: 1
      const double Gip0 = W27[ip]*G0;
      
      // Number of operations: 1
      const double Gip1 = W27[ip]*G1;
      
      // Number of operations: 1
      const double Gip2 = W27[ip]*G2;
      
      // Number of operations: 1
      const double Gip3 = W27[ip]*G3;
      
      // Number of operations: 5
      const double Gip4 = W27[ip]*(G4 + F0*(G5 + F0*G6));
      
      // Number of operations: 1
      const double Gip5 = W27[ip]*G7;
      
      // Number of operations: 1
      const double Gip6 = W27[ip]*G8;
      
      // Number of operations: 1
      const double Gip7 = W27[ip]*G9;
      
      // Number of operations: 1
      const double Gip8 = W27[ip]*det;
      
      // Number of operations: 1
      const double Gip9 = W27[ip]*G10;
      
      // Number of operations: 1
      const double Gip10 = W27[ip]*G11;
      
      // Number of operations: 1
      const double Gip11 = W27[ip]*G12;
      
      // Number of operations: 1
      const double Gip12 = W27[ip]*G13;
      
      // Number of operations: 1
      const double Gip13 = W27[ip]*G14;
      
      
      // Number of operations for primary indices = 216
      for (unsigned int j = 0; j < 2; j++)
      {
        for (unsigned int k = 0; k < 2; k++)
        {
          // Number of operations to compute entry = 3
          A[nzc1[j]*8 + nzc1[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip0;
          // Number of operations to compute entry = 3
          A[nzc4[j]*8 + nzc3[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip1;
          // Number of operations to compute entry = 3
          A[nzc0[j]*8 + nzc2[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip2;
          // Number of operations to compute entry = 3
          A[nzc1[j]*8 + nzc2[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip3;
          // Number of operations to compute entry = 3
          A[nzc2[j]*8 + nzc1[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip3;
          // Number of operations to compute entry = 3
          A[nzc3[j]*8 + nzc5[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip5;
          // Number of operations to compute entry = 3
          A[nzc4[j]*8 + nzc4[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip6;
          // Number of operations to compute entry = 3
          A[nzc5[j]*8 + nzc4[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip7;
          // Number of operations to compute entry = 3
          A[nzc4[j]*8 + nzc5[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip7;
          // Number of operations to compute entry = 3
          A[nzc2[j]*8 + nzc2[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip9;
          // Number of operations to compute entry = 3
          A[nzc5[j]*8 + nzc3[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip5;
          // Number of operations to compute entry = 3
          A[nzc1[j]*8 + nzc0[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip10;
          // Number of operations to compute entry = 3
          A[nzc5[j]*8 + nzc5[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip11;
          // Number of operations to compute entry = 3
          A[nzc3[j]*8 + nzc4[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip1;
          // Number of operations to compute entry = 3
          A[nzc0[j]*8 + nzc0[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip12;
          // Number of operations to compute entry = 3
          A[nzc3[j]*8 + nzc3[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip13;
          // Number of operations to compute entry = 3
          A[nzc2[j]*8 + nzc0[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip2;
          // Number of operations to compute entry = 3
          A[nzc0[j]*8 + nzc1[k]] += FE0_C1_D001[ip][j]*FE0_C1_D001[ip][k]*Gip10;
        }// end loop over 'k'
      }// end loop over 'j'
      
      // Number of operations for primary indices = 144
      for (unsigned int j = 0; j < 4; j++)
      {
        for (unsigned int k = 0; k < 4; k++)
        {
          // Number of operations to compute entry = 3
          A[nzc6[j]*8 + nzc6[k]] += FE0_C1[ip][j]*FE0_C1[ip][k]*Gip4;
          // Number of operations to compute entry = 3
          A[nzc7[j]*8 + nzc6[k]] += FE0_C1[ip][j]*FE0_C1[ip][k]*Gip8;
          // Number of operations to compute entry = 3
          A[nzc6[j]*8 + nzc7[k]] += FE0_C1[ip][j]*FE0_C1[ip][k]*Gip8;
        }// end loop over 'k'
      }// end loop over 'j'
    }// end loop over 'ip'
}

/// Constructor
cahnhilliard3d_0_cell_integral_0::cahnhilliard3d_0_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_0_cell_integral_0::~cahnhilliard3d_0_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void cahnhilliard3d_0_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Reset values of the element tensor block
    A[0] = 0;
    A[1] = 0;
    A[2] = 0;
    A[3] = 0;
    A[4] = 0;
    A[5] = 0;
    A[6] = 0;
    A[7] = 0;
    A[8] = 0;
    A[9] = 0;
    A[10] = 0;
    A[11] = 0;
    A[12] = 0;
    A[13] = 0;
    A[14] = 0;
    A[15] = 0;
    A[16] = 0;
    A[17] = 0;
    A[18] = 0;
    A[19] = 0;
    A[20] = 0;
    A[21] = 0;
    A[22] = 0;
    A[23] = 0;
    A[24] = 0;
    A[25] = 0;
    A[26] = 0;
    A[27] = 0;
    A[28] = 0;
    A[29] = 0;
    A[30] = 0;
    A[31] = 0;
    A[32] = 0;
    A[33] = 0;
    A[34] = 0;
    A[35] = 0;
    A[36] = 0;
    A[37] = 0;
    A[38] = 0;
    A[39] = 0;
    A[40] = 0;
    A[41] = 0;
    A[42] = 0;
    A[43] = 0;
    A[44] = 0;
    A[45] = 0;
    A[46] = 0;
    A[47] = 0;
    A[48] = 0;
    A[49] = 0;
    A[50] = 0;
    A[51] = 0;
    A[52] = 0;
    A[53] = 0;
    A[54] = 0;
    A[55] = 0;
    A[56] = 0;
    A[57] = 0;
    A[58] = 0;
    A[59] = 0;
    A[60] = 0;
    A[61] = 0;
    A[62] = 0;
    A[63] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c);
}

/// Constructor
cahnhilliard3d_form_0::cahnhilliard3d_form_0() : ufc::form()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_form_0::~cahnhilliard3d_form_0()
{
    // Do nothing
}

/// Return a string identifying the form
const char* cahnhilliard3d_form_0::signature() const
{
    return "Form([Integral(Sum(Sum(Product(Constant(Cell('tetrahedron', 1), 3), IndexSum(Product(Indexed(ComponentTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((Index(0),), {Index(0): 3})), MultiIndex((FixedIndex(0),), {})), MultiIndex((Index(0),), {Index(0): 3})), MultiIndex((Index(1),), {Index(1): 3})), Indexed(ComponentTensor(Product(Constant(Cell('tetrahedron', 1), 4), Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((Index(2),), {Index(2): 3})), MultiIndex((FixedIndex(0),), {}))), MultiIndex((Index(2),), {Index(2): 3})), MultiIndex((Index(1),), {Index(1): 3}))), MultiIndex((Index(1),), {Index(1): 3}))), Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 2})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))), Sum(Product(IntValue(-1, (), (), {}), Product(Constant(Cell('tetrahedron', 1), 1), IndexSum(Product(Indexed(ComponentTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((Index(3),), {Index(3): 3})), MultiIndex((FixedIndex(1),), {})), MultiIndex((Index(3),), {Index(3): 3})), MultiIndex((Index(4),), {Index(4): 3})), Indexed(ComponentTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((Index(5),), {Index(5): 3})), MultiIndex((FixedIndex(1),), {})), MultiIndex((Index(5),), {Index(5): 3})), MultiIndex((Index(4),), {Index(4): 3}))), MultiIndex((Index(4),), {Index(4): 3})))), Sum(Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((FixedIndex(0),), {FixedIndex(0): 2}))), Product(IntValue(-1, (), (), {}), Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Product(Constant(Cell('tetrahedron', 1), 2), Sum(Product(Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Product(FloatValue(2.0, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))), Sum(Product(IntValue(-1, (), (), {}), Product(Product(FloatValue(8.0, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2}))), Sum(FloatValue(1.0, (), (), {}), Product(IntValue(-1, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))))), Product(Product(FloatValue(2.0, (), (), {}), Sum(FloatValue(1.0, (), (), {}), Product(IntValue(-1, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2}))))), Sum(FloatValue(1.0, (), (), {}), Product(IntValue(-1, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))))))))))))), Measure('cell', 0, None))])";
}

/// Return the rank of the global tensor (r)
unsigned int cahnhilliard3d_form_0::rank() const
{
    return 2;
}

/// Return the number of coefficients (n)
unsigned int cahnhilliard3d_form_0::num_coefficients() const
{
    return 5;
}

/// Return the number of cell integrals
unsigned int cahnhilliard3d_form_0::num_cell_integrals() const
{
    return 1;
}
  
/// Return the number of exterior facet integrals
unsigned int cahnhilliard3d_form_0::num_exterior_facet_integrals() const
{
    return 0;
}
  
/// Return the number of interior facet integrals
unsigned int cahnhilliard3d_form_0::num_interior_facet_integrals() const
{
    return 0;
}

/// Create a new finite element for argument function i
ufc::finite_element* cahnhilliard3d_form_0::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_finite_element_0();
      break;
    case 1:
      return new cahnhilliard3d_0_finite_element_1();
      break;
    case 2:
      return new cahnhilliard3d_0_finite_element_2();
      break;
    case 3:
      return new cahnhilliard3d_0_finite_element_3();
      break;
    case 4:
      return new cahnhilliard3d_0_finite_element_4();
      break;
    case 5:
      return new cahnhilliard3d_0_finite_element_5();
      break;
    case 6:
      return new cahnhilliard3d_0_finite_element_6();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* cahnhilliard3d_form_0::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_0_dof_map_0();
      break;
    case 1:
      return new cahnhilliard3d_0_dof_map_1();
      break;
    case 2:
      return new cahnhilliard3d_0_dof_map_2();
      break;
    case 3:
      return new cahnhilliard3d_0_dof_map_3();
      break;
    case 4:
      return new cahnhilliard3d_0_dof_map_4();
      break;
    case 5:
      return new cahnhilliard3d_0_dof_map_5();
      break;
    case 6:
      return new cahnhilliard3d_0_dof_map_6();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* cahnhilliard3d_form_0::create_cell_integral(unsigned int i) const
{
    return new cahnhilliard3d_0_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* cahnhilliard3d_form_0::create_exterior_facet_integral(unsigned int i) const
{
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* cahnhilliard3d_form_0::create_interior_facet_integral(unsigned int i) const
{
    return 0;
}


/// Constructor
cahnhilliard3d_1_finite_element_0_0::cahnhilliard3d_1_finite_element_0_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_0_0::~cahnhilliard3d_1_finite_element_0_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_0_0::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_0_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_0_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_0_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_0_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_0_0::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_1_finite_element_0_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_0_0::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_0_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_0_0::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_1_finite_element_0_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_0_0::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_0_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_0_0::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_0_0();
}


/// Constructor
cahnhilliard3d_1_finite_element_0_1::cahnhilliard3d_1_finite_element_0_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_0_1::~cahnhilliard3d_1_finite_element_0_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_0_1::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_0_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_0_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_0_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_0_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_0_1::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_1_finite_element_0_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_0_1::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_0_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_0_1::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_1_finite_element_0_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_0_1::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_0_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_0_1::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_0_1();
}


/// Constructor
cahnhilliard3d_1_finite_element_0::cahnhilliard3d_1_finite_element_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_0::~cahnhilliard3d_1_finite_element_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_0::signature() const
{
    return "MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_0::space_dimension() const
{
    return 8;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_0::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_0::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_0::evaluate_basis(unsigned int i,
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
    values[0] = 0;
    values[1] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_0::evaluate_basis_derivatives(unsigned int i,
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
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
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_0::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[8][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[8][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[8][1][2] = {{{1, 0}}, {{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}, {{0, 1}}, {{0, 1}}, {{0, 1}}};
    
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
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void cahnhilliard3d_1_finite_element_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_0::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[2] = dof_values[1];
    vertex_values[4] = dof_values[2];
    vertex_values[6] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[3] = dof_values[5];
    vertex_values[5] = dof_values[6];
    vertex_values[7] = dof_values[7];
}

/// Return the number of sub elements (for a mixed element)
unsigned int cahnhilliard3d_1_finite_element_0::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_0::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_finite_element_0_0();
      break;
    case 1:
      return new cahnhilliard3d_1_finite_element_0_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_1_finite_element_1_0::cahnhilliard3d_1_finite_element_1_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_1_0::~cahnhilliard3d_1_finite_element_1_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_1_0::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_1_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_1_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_1_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_1_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_1_0::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_1_finite_element_1_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_1_0::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_1_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_1_0::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_1_finite_element_1_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_1_0::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_1_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_1_0::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_1_0();
}


/// Constructor
cahnhilliard3d_1_finite_element_1_1::cahnhilliard3d_1_finite_element_1_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_1_1::~cahnhilliard3d_1_finite_element_1_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_1_1::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_1_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_1_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_1_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_1_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_1_1::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_1_finite_element_1_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_1_1::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_1_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_1_1::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_1_finite_element_1_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_1_1::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_1_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_1_1::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_1_1();
}


/// Constructor
cahnhilliard3d_1_finite_element_1::cahnhilliard3d_1_finite_element_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_1::~cahnhilliard3d_1_finite_element_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_1::signature() const
{
    return "MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_1::space_dimension() const
{
    return 8;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_1::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_1::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_1::evaluate_basis(unsigned int i,
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
    values[0] = 0;
    values[1] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_1::evaluate_basis_derivatives(unsigned int i,
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
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
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_1::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[8][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[8][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[8][1][2] = {{{1, 0}}, {{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}, {{0, 1}}, {{0, 1}}, {{0, 1}}};
    
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
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void cahnhilliard3d_1_finite_element_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_1::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[2] = dof_values[1];
    vertex_values[4] = dof_values[2];
    vertex_values[6] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[3] = dof_values[5];
    vertex_values[5] = dof_values[6];
    vertex_values[7] = dof_values[7];
}

/// Return the number of sub elements (for a mixed element)
unsigned int cahnhilliard3d_1_finite_element_1::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_1::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_finite_element_1_0();
      break;
    case 1:
      return new cahnhilliard3d_1_finite_element_1_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_1_finite_element_2_0::cahnhilliard3d_1_finite_element_2_0() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_2_0::~cahnhilliard3d_1_finite_element_2_0()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_2_0::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_2_0::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_2_0::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_2_0::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_2_0::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_2_0::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_1_finite_element_2_0::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_2_0::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_2_0::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_2_0::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_1_finite_element_2_0::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_2_0::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_2_0::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_2_0::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_2_0();
}


/// Constructor
cahnhilliard3d_1_finite_element_2_1::cahnhilliard3d_1_finite_element_2_1() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_2_1::~cahnhilliard3d_1_finite_element_2_1()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_2_1::signature() const
{
    return "FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_2_1::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_2_1::space_dimension() const
{
    return 4;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_2_1::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_2_1::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_2_1::evaluate_basis(unsigned int i,
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
void cahnhilliard3d_1_finite_element_2_1::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_2_1::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_2_1::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_2_1::evaluate_dof(unsigned int i,
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
void cahnhilliard3d_1_finite_element_2_1::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_2_1::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_2_1::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_2_1::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_2_1();
}


/// Constructor
cahnhilliard3d_1_finite_element_2::cahnhilliard3d_1_finite_element_2() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_2::~cahnhilliard3d_1_finite_element_2()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_2::signature() const
{
    return "MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_2::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_2::space_dimension() const
{
    return 8;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_2::value_rank() const
{
    return 1;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_2::value_dimension(unsigned int i) const
{
    return 2;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_2::evaluate_basis(unsigned int i,
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
    values[0] = 0;
    values[1] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[0] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Extract relevant coefficients
      const double coeff0_0 =   coefficients0[dof][0];
      const double coeff0_1 =   coefficients0[dof][1];
      const double coeff0_2 =   coefficients0[dof][2];
      const double coeff0_3 =   coefficients0[dof][3];
    
      // Compute value(s)
      values[1] = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3;
    }
    
}

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_2::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_2::evaluate_basis_derivatives(unsigned int i,
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    if (0 <= i && i <= 3)
    {
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
    
    if (4 <= i && i <= 7)
    {
      // Map degree of freedom to element degree of freedom
      const unsigned int dof = i - 4;
    
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
      const static double coefficients0[4][4] =   \
      {{0.288675134594813, -0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0.182574185835055, -0.105409255338946, -0.074535599249993},
      {0.288675134594813, 0, 0.210818510677892, -0.074535599249993},
      {0.288675134594813, 0, 0, 0.223606797749979}};
    
      // Interesting (new) part
      // Tables of derivatives of the polynomial base (transpose)
      const static double dmats0[4][4] =   \
      {{0, 0, 0, 0},
      {6.32455532033676, 0, 0, 0},
      {0, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats1[4][4] =   \
      {{0, 0, 0, 0},
      {3.16227766016838, 0, 0, 0},
      {5.47722557505166, 0, 0, 0},
      {0, 0, 0, 0}};
    
      const static double dmats2[4][4] =   \
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
          values[num_derivatives + row] += transform[row][col]*derivatives[col];
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
    
}

/// Evaluate order n derivatives of all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_2::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_2::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[8][1][3] = {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}, {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    const static double W[8][1] = {{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[8][1][2] = {{{1, 0}}, {{1, 0}}, {{1, 0}}, {{1, 0}}, {{0, 1}}, {{0, 1}}, {{0, 1}}, {{0, 1}}};
    
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
    double values[2];
    f.evaluate(values, y, c);
    
    // Map function values using appropriate mapping
    // Affine map: Do nothing
    
    // Note that we do not map the weights (yet).
    
    // Take directional components
    for(int k = 0; k < 2; k++)
      result += values[k]*D[i][0][k];
    // Multiply by weights 
    result *= W[i][0];
    
    return result;
}

/// Evaluate linear functionals for all dofs on the function f
void cahnhilliard3d_1_finite_element_2::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_2::interpolate_vertex_values(double* vertex_values,
                                              const double* dof_values,
                                              const ufc::cell& c) const
{
    // Evaluate at vertices and use affine mapping
    vertex_values[0] = dof_values[0];
    vertex_values[2] = dof_values[1];
    vertex_values[4] = dof_values[2];
    vertex_values[6] = dof_values[3];
    // Evaluate at vertices and use affine mapping
    vertex_values[1] = dof_values[4];
    vertex_values[3] = dof_values[5];
    vertex_values[5] = dof_values[6];
    vertex_values[7] = dof_values[7];
}

/// Return the number of sub elements (for a mixed element)
unsigned int cahnhilliard3d_1_finite_element_2::num_sub_elements() const
{
    return 2;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_2::create_sub_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_finite_element_2_0();
      break;
    case 1:
      return new cahnhilliard3d_1_finite_element_2_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_1_finite_element_3::cahnhilliard3d_1_finite_element_3() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_3::~cahnhilliard3d_1_finite_element_3()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_3::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_3::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_3::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_3::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_3::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_3::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_3::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_3::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_3::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_3::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_1_finite_element_3::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_3::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_3::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_3::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_3();
}


/// Constructor
cahnhilliard3d_1_finite_element_4::cahnhilliard3d_1_finite_element_4() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_4::~cahnhilliard3d_1_finite_element_4()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_4::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_4::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_4::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_4::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_4::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_4::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_4::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_4::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_4::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_4::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_1_finite_element_4::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_4::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_4::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_4::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_4();
}


/// Constructor
cahnhilliard3d_1_finite_element_5::cahnhilliard3d_1_finite_element_5() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_5::~cahnhilliard3d_1_finite_element_5()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_5::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_5::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_5::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_5::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_5::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_5::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_5::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_5::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_5::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_5::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_1_finite_element_5::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_5::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_5::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_5::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_5();
}


/// Constructor
cahnhilliard3d_1_finite_element_6::cahnhilliard3d_1_finite_element_6() : ufc::finite_element()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_finite_element_6::~cahnhilliard3d_1_finite_element_6()
{
    // Do nothing
}

/// Return a string identifying the finite element
const char* cahnhilliard3d_1_finite_element_6::signature() const
{
    return "FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return the cell shape
ufc::shape cahnhilliard3d_1_finite_element_6::cell_shape() const
{
    return ufc::tetrahedron;
}

/// Return the dimension of the finite element function space
unsigned int cahnhilliard3d_1_finite_element_6::space_dimension() const
{
    return 1;
}

/// Return the rank of the value space
unsigned int cahnhilliard3d_1_finite_element_6::value_rank() const
{
    return 0;
}

/// Return the dimension of the value space for axis i
unsigned int cahnhilliard3d_1_finite_element_6::value_dimension(unsigned int i) const
{
    return 1;
}

/// Evaluate basis function i at given point in cell
void cahnhilliard3d_1_finite_element_6::evaluate_basis(unsigned int i,
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

/// Evaluate all basis functions at given point in cell
void cahnhilliard3d_1_finite_element_6::evaluate_basis_all(double* values,
                                       const double* coordinates,
                                       const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis() is not yet implemented.");
}

/// Evaluate order n derivatives of basis function i at given point in cell
void cahnhilliard3d_1_finite_element_6::evaluate_basis_derivatives(unsigned int i,
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
void cahnhilliard3d_1_finite_element_6::evaluate_basis_derivatives_all(unsigned int n,
                                                   double* values,
                                                   const double* coordinates,
                                                   const ufc::cell& c) const
{
    throw std::runtime_error("The vectorised version of evaluate_basis_derivatives() is not yet implemented.");
}

/// Evaluate linear functional for dof i on the function f
double cahnhilliard3d_1_finite_element_6::evaluate_dof(unsigned int i,
                                   const ufc::function& f,
                                   const ufc::cell& c) const
{
    // The reference points, direction and weights:
    const static double X[1][1][3] = {{{0.25, 0.25, 0.25}}};
    const static double W[1][1] = {{1}};
    const static double D[1][1][1] = {{{1}}};
    
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
void cahnhilliard3d_1_finite_element_6::evaluate_dofs(double* values,
                                  const ufc::function& f,
                                  const ufc::cell& c) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Interpolate vertex values from dof values
void cahnhilliard3d_1_finite_element_6::interpolate_vertex_values(double* vertex_values,
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
unsigned int cahnhilliard3d_1_finite_element_6::num_sub_elements() const
{
    return 1;
}

/// Create a new finite element for sub element i (for a mixed element)
ufc::finite_element* cahnhilliard3d_1_finite_element_6::create_sub_element(unsigned int i) const
{
    return new cahnhilliard3d_1_finite_element_6();
}

/// Constructor
cahnhilliard3d_1_dof_map_0_0::cahnhilliard3d_1_dof_map_0_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_0_0::~cahnhilliard3d_1_dof_map_0_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_0_0::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_0_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_0_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_0_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_0_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_0_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_0_0::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_0_0::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_0_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_0_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_0_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_0_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_0_0::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_0_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_0_0::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_1_dof_map_0_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_0_0::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_0_0();
}


/// Constructor
cahnhilliard3d_1_dof_map_0_1::cahnhilliard3d_1_dof_map_0_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_0_1::~cahnhilliard3d_1_dof_map_0_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_0_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_0_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_0_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_0_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_0_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_0_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_0_1::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_0_1::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_0_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_0_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_0_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_0_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_0_1::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_0_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_0_1::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_1_dof_map_0_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_0_1::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_0_1();
}


/// Constructor
cahnhilliard3d_1_dof_map_0::cahnhilliard3d_1_dof_map_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_0::~cahnhilliard3d_1_dof_map_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_0::signature() const
{
    return "FFC dof map for MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_0::local_dimension(const ufc::cell& c) const
{
    return 8;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_0::max_local_dimension() const
{
    return 8;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_0::num_facet_dofs() const
{
    return 6;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_0::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_0::tabulate_coordinates(double** coordinates,
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
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_0::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_0::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_dof_map_0_0();
      break;
    case 1:
      return new cahnhilliard3d_1_dof_map_0_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_1_dof_map_1_0::cahnhilliard3d_1_dof_map_1_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_1_0::~cahnhilliard3d_1_dof_map_1_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_1_0::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_1_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_1_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_1_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_1_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_1_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_1_0::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_1_0::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_1_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_1_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_1_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_1_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_1_0::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_1_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_1_0::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_1_dof_map_1_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_1_0::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_1_0();
}


/// Constructor
cahnhilliard3d_1_dof_map_1_1::cahnhilliard3d_1_dof_map_1_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_1_1::~cahnhilliard3d_1_dof_map_1_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_1_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_1_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_1_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_1_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_1_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_1_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_1_1::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_1_1::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_1_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_1_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_1_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_1_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_1_1::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_1_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_1_1::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_1_dof_map_1_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_1_1::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_1_1();
}


/// Constructor
cahnhilliard3d_1_dof_map_1::cahnhilliard3d_1_dof_map_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_1::~cahnhilliard3d_1_dof_map_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_1::signature() const
{
    return "FFC dof map for MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_1::local_dimension(const ufc::cell& c) const
{
    return 8;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_1::max_local_dimension() const
{
    return 8;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_1::num_facet_dofs() const
{
    return 6;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_1::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_1::tabulate_coordinates(double** coordinates,
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
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_1::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_1::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_dof_map_1_0();
      break;
    case 1:
      return new cahnhilliard3d_1_dof_map_1_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_1_dof_map_2_0::cahnhilliard3d_1_dof_map_2_0() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_2_0::~cahnhilliard3d_1_dof_map_2_0()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_2_0::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_2_0::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_2_0::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_2_0::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_2_0::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_2_0::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_2_0::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_2_0::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_2_0::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_2_0::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_2_0::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_2_0::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_2_0::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_2_0::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_2_0::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_1_dof_map_2_0::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_2_0::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_2_0();
}


/// Constructor
cahnhilliard3d_1_dof_map_2_1::cahnhilliard3d_1_dof_map_2_1() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_2_1::~cahnhilliard3d_1_dof_map_2_1()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_2_1::signature() const
{
    return "FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_2_1::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_2_1::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_2_1::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_2_1::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_2_1::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_2_1::local_dimension(const ufc::cell& c) const
{
    return 4;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_2_1::max_local_dimension() const
{
    return 4;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_2_1::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_2_1::num_facet_dofs() const
{
    return 3;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_2_1::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_2_1::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_2_1::tabulate_facet_dofs(unsigned int* dofs,
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

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_2_1::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_2_1::tabulate_coordinates(double** coordinates,
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
unsigned int cahnhilliard3d_1_dof_map_2_1::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_2_1::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_2_1();
}


/// Constructor
cahnhilliard3d_1_dof_map_2::cahnhilliard3d_1_dof_map_2() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_2::~cahnhilliard3d_1_dof_map_2()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_2::signature() const
{
    return "FFC dof map for MixedElement([FiniteElement('Lagrange', 'tetrahedron', 1), FiniteElement('Lagrange', 'tetrahedron', 1)])";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_2::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_2::init_mesh(const ufc::mesh& m)
{
    __global_dimension = 2*m.num_entities[0];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_2::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_2::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_2::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_2::local_dimension(const ufc::cell& c) const
{
    return 8;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_2::max_local_dimension() const
{
    return 8;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_2::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_2::num_facet_dofs() const
{
    return 6;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_2::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_2::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
    dofs[3] = c.entity_indices[0][3];
    unsigned int offset = m.num_entities[0];
    dofs[4] = offset + c.entity_indices[0][0];
    dofs[5] = offset + c.entity_indices[0][1];
    dofs[6] = offset + c.entity_indices[0][2];
    dofs[7] = offset + c.entity_indices[0][3];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_2::tabulate_facet_dofs(unsigned int* dofs,
                                        unsigned int facet) const
{
    switch ( facet )
    {
    case 0:
      dofs[0] = 1;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 5;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 6;
      dofs[5] = 7;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 3;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 7;
      break;
    case 3:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 2;
      dofs[3] = 4;
      dofs[4] = 5;
      dofs[5] = 6;
      break;
    }
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void cahnhilliard3d_1_dof_map_2::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_2::tabulate_coordinates(double** coordinates,
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
    coordinates[4][0] = x[0][0];
    coordinates[4][1] = x[0][1];
    coordinates[4][2] = x[0][2];
    coordinates[5][0] = x[1][0];
    coordinates[5][1] = x[1][1];
    coordinates[5][2] = x[1][2];
    coordinates[6][0] = x[2][0];
    coordinates[6][1] = x[2][1];
    coordinates[6][2] = x[2][2];
    coordinates[7][0] = x[3][0];
    coordinates[7][1] = x[3][1];
    coordinates[7][2] = x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_2::num_sub_dof_maps() const
{
    return 2;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_2::create_sub_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_dof_map_2_0();
      break;
    case 1:
      return new cahnhilliard3d_1_dof_map_2_1();
      break;
    }
    return 0;
}


/// Constructor
cahnhilliard3d_1_dof_map_3::cahnhilliard3d_1_dof_map_3() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_3::~cahnhilliard3d_1_dof_map_3()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_3::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_3::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_3::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_3::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_3::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_3::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_3::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_3::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_3::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_3::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_3::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_3::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_3::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_1_dof_map_3::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_3::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_3::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_3::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_3();
}


/// Constructor
cahnhilliard3d_1_dof_map_4::cahnhilliard3d_1_dof_map_4() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_4::~cahnhilliard3d_1_dof_map_4()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_4::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_4::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_4::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_4::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_4::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_4::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_4::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_4::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_4::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_4::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_4::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_4::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_4::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_1_dof_map_4::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_4::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_4::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_4::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_4();
}


/// Constructor
cahnhilliard3d_1_dof_map_5::cahnhilliard3d_1_dof_map_5() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_5::~cahnhilliard3d_1_dof_map_5()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_5::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_5::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_5::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_5::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_5::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_5::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_5::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_5::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_5::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_5::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_5::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_5::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_5::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_1_dof_map_5::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_5::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_5::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_5::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_5();
}


/// Constructor
cahnhilliard3d_1_dof_map_6::cahnhilliard3d_1_dof_map_6() : ufc::dof_map()
{
    __global_dimension = 0;
}

/// Destructor
cahnhilliard3d_1_dof_map_6::~cahnhilliard3d_1_dof_map_6()
{
    // Do nothing
}

/// Return a string identifying the dof map
const char* cahnhilliard3d_1_dof_map_6::signature() const
{
    return "FFC dof map for FiniteElement('Discontinuous Lagrange', 'tetrahedron', 0)";
}

/// Return true iff mesh entities of topological dimension d are needed
bool cahnhilliard3d_1_dof_map_6::needs_mesh_entities(unsigned int d) const
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
bool cahnhilliard3d_1_dof_map_6::init_mesh(const ufc::mesh& m)
{
    __global_dimension = m.num_entities[3];
    return false;
}

/// Initialize dof map for given cell
void cahnhilliard3d_1_dof_map_6::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
    // Do nothing
}

/// Finish initialization of dof map for cells
void cahnhilliard3d_1_dof_map_6::init_cell_finalize()
{
    // Do nothing
}

/// Return the dimension of the global finite element function space
unsigned int cahnhilliard3d_1_dof_map_6::global_dimension() const
{
    return __global_dimension;
}

/// Return the dimension of the local finite element function space for a cell
unsigned int cahnhilliard3d_1_dof_map_6::local_dimension(const ufc::cell& c) const
{
    return 1;
}

/// Return the maximum dimension of the local finite element function space
unsigned int cahnhilliard3d_1_dof_map_6::max_local_dimension() const
{
    return 1;
}

// Return the geometric dimension of the coordinates this dof map provides
unsigned int cahnhilliard3d_1_dof_map_6::geometric_dimension() const
{
    return 3;
}

/// Return the number of dofs on each cell facet
unsigned int cahnhilliard3d_1_dof_map_6::num_facet_dofs() const
{
    return 0;
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int cahnhilliard3d_1_dof_map_6::num_entity_dofs(unsigned int d) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the local-to-global mapping of dofs on a cell
void cahnhilliard3d_1_dof_map_6::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
    dofs[0] = c.entity_indices[3][0];
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void cahnhilliard3d_1_dof_map_6::tabulate_facet_dofs(unsigned int* dofs,
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
void cahnhilliard3d_1_dof_map_6::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
    throw std::runtime_error("Not implemented (introduced in UFC v1.1).");
}

/// Tabulate the coordinates of all dofs on a cell
void cahnhilliard3d_1_dof_map_6::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
    const double * const * x = c.coordinates;
    coordinates[0][0] = 0.25*x[0][0] + 0.25*x[1][0] + 0.25*x[2][0] + 0.25*x[3][0];
    coordinates[0][1] = 0.25*x[0][1] + 0.25*x[1][1] + 0.25*x[2][1] + 0.25*x[3][1];
    coordinates[0][2] = 0.25*x[0][2] + 0.25*x[1][2] + 0.25*x[2][2] + 0.25*x[3][2];
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int cahnhilliard3d_1_dof_map_6::num_sub_dof_maps() const
{
    return 1;
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* cahnhilliard3d_1_dof_map_6::create_sub_dof_map(unsigned int i) const
{
    return new cahnhilliard3d_1_dof_map_6();
}


/// Constructor
cahnhilliard3d_1_cell_integral_0_quadrature::cahnhilliard3d_1_cell_integral_0_quadrature() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_cell_integral_0_quadrature::~cahnhilliard3d_1_cell_integral_0_quadrature()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void cahnhilliard3d_1_cell_integral_0_quadrature::tabulate_tensor(double* A,
                                    const double * const * w,
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
    const double Jinv_00 = d_00 / detJ;
    const double Jinv_01 = d_10 / detJ;
    const double Jinv_02 = d_20 / detJ;
    const double Jinv_10 = d_01 / detJ;
    const double Jinv_11 = d_11 / detJ;
    const double Jinv_12 = d_21 / detJ;
    const double Jinv_20 = d_02 / detJ;
    const double Jinv_21 = d_12 / detJ;
    const double Jinv_22 = d_22 / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    
    // Array of quadrature weights
    const static double W27[27] = {0.00877047492965105, 0.00816265076654668, 0.0016716811314837, 0.0100061425721761, 0.00931268237947044, 0.00190720341498178, 0.00304787709051818, 0.00283664869563092, 0.000580935315837384, 0.0140327598874417, 0.0130602412264747, 0.00267468981037392, 0.0160098281154818, 0.0149002918071527, 0.00305152546397085, 0.0048766033448291, 0.00453863791300947, 0.000929496505339815, 0.00877047492965105, 0.00816265076654668, 0.0016716811314837, 0.0100061425721761, 0.00931268237947044, 0.00190720341498178, 0.00304787709051818, 0.00283664869563092, 0.000580935315837384};
    
    
    const static double FE0_C1_D001[27][2] = \
    {{-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1},
    {-1, 1}};
    // Array of non-zero columns
    static const unsigned int nzc0[2] = {4, 7};
    // Array of non-zero columns
    static const unsigned int nzc1[2] = {4, 6};
    // Array of non-zero columns
    static const unsigned int nzc2[2] = {4, 5};
    // Array of non-zero columns
    static const unsigned int nzc3[2] = {0, 1};
    // Array of non-zero columns
    static const unsigned int nzc4[2] = {0, 2};
    // Array of non-zero columns
    static const unsigned int nzc5[2] = {0, 3};
    
    const static double FE0_C1[27][4] = \
    {{0.749664528221693, 0.0952198798417149, 0.0821215678634426, 0.0729940240731498},
    {0.528074388273447, 0.0670742417520586, 0.0578476039361428, 0.347003766038352},
    {0.238563056650491, 0.0303014811742758, 0.0261332522867349, 0.705002209888498},
    {0.485731727037113, 0.0616960186091465, 0.379578230280591, 0.0729940240731498},
    {0.342156357895961, 0.0434595556538024, 0.267380320411885, 0.347003766038352},
    {0.154572667042115, 0.0196333029354846, 0.120791820133903, 0.705002209888498},
    {0.174656645238399, 0.0221843026408197, 0.730165028047632, 0.0729940240731498},
    {0.123030632529655, 0.0156269392579017, 0.514338662174092, 0.347003766038352},
    {0.0555803583920822, 0.00705963113955482, 0.232357800579865, 0.705002209888498},
    {0.422442204031704, 0.422442204031704, 0.0821215678634426, 0.0729940240731498},
    {0.297574315012753, 0.297574315012753, 0.0578476039361428, 0.347003766038352},
    {0.134432268912383, 0.134432268912383, 0.0261332522867349, 0.705002209888498},
    {0.27371387282313, 0.27371387282313, 0.379578230280591, 0.0729940240731498},
    {0.192807956774882, 0.192807956774882, 0.267380320411885, 0.347003766038352},
    {0.0871029849887995, 0.0871029849887996, 0.120791820133903, 0.705002209888498},
    {0.0984204739396093, 0.0984204739396094, 0.730165028047632, 0.0729940240731498},
    {0.0693287858937781, 0.0693287858937782, 0.514338662174092, 0.347003766038352},
    {0.0313199947658185, 0.0313199947658185, 0.232357800579865, 0.705002209888498},
    {0.0952198798417148, 0.749664528221693, 0.0821215678634426, 0.0729940240731498},
    {0.0670742417520584, 0.528074388273447, 0.0578476039361428, 0.347003766038352},
    {0.0303014811742757, 0.238563056650491, 0.0261332522867349, 0.705002209888498},
    {0.0616960186091464, 0.485731727037113, 0.379578230280591, 0.0729940240731498},
    {0.0434595556538024, 0.342156357895961, 0.267380320411885, 0.347003766038352},
    {0.0196333029354845, 0.154572667042115, 0.120791820133903, 0.705002209888498},
    {0.0221843026408197, 0.174656645238399, 0.730165028047632, 0.0729940240731498},
    {0.0156269392579016, 0.123030632529655, 0.514338662174092, 0.347003766038352},
    {0.00705963113955483, 0.0555803583920822, 0.232357800579865, 0.705002209888498}};
    // Array of non-zero columns
    static const unsigned int nzc6[4] = {4, 5, 6, 7};
    // Array of non-zero columns
    static const unsigned int nzc7[4] = {0, 1, 2, 3};
    
    // Number of operations to compute geometry constants: 189
    const double G0 =  - det*w[2][0]*(Jinv_10*Jinv_20 + Jinv_11*Jinv_21 + Jinv_12*Jinv_22);
    const double G1 =  - det*w[2][0]*(Jinv_10*Jinv_10 + Jinv_11*Jinv_11 + Jinv_12*Jinv_12);
    const double G2 =  - det*w[2][0]*(Jinv_00*Jinv_10 + Jinv_01*Jinv_11 + Jinv_02*Jinv_12);
    const double G3 = det*w[4][0]*(Jinv_10*Jinv_20 + Jinv_11*Jinv_21 + Jinv_12*Jinv_22 - w[5][0]*(Jinv_10*Jinv_20 + Jinv_11*Jinv_21 + Jinv_12*Jinv_22));
    const double G4 = det*w[4][0]*(Jinv_10*Jinv_10 + Jinv_11*Jinv_11 + Jinv_12*Jinv_12 - w[5][0]*(Jinv_10*Jinv_10 + Jinv_11*Jinv_11 + Jinv_12*Jinv_12));
    const double G5 = det*w[4][0]*w[5][0]*(Jinv_10*Jinv_10 + Jinv_11*Jinv_11 + Jinv_12*Jinv_12);
    const double G6 = det*w[4][0]*(Jinv_00*Jinv_10 + Jinv_01*Jinv_11 + Jinv_02*Jinv_12 - w[5][0]*(Jinv_00*Jinv_10 + Jinv_01*Jinv_11 + Jinv_02*Jinv_12));
    const double G7 = det*w[4][0]*w[5][0]*(Jinv_10*Jinv_20 + Jinv_11*Jinv_21 + Jinv_12*Jinv_22);
    const double G8 = det*w[4][0]*w[5][0]*(Jinv_00*Jinv_10 + Jinv_01*Jinv_11 + Jinv_02*Jinv_12);
    const double G9 =  - det*w[2][0]*(Jinv_00*Jinv_20 + Jinv_01*Jinv_21 + Jinv_02*Jinv_22);
    const double G10 =  - det*w[2][0]*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01 + Jinv_02*Jinv_02);
    const double G11 = det*w[4][0]*(Jinv_20*Jinv_20 + Jinv_21*Jinv_21 + Jinv_22*Jinv_22 - w[5][0]*(Jinv_20*Jinv_20 + Jinv_21*Jinv_21 + Jinv_22*Jinv_22));
    const double G12 = det*w[4][0]*(Jinv_00*Jinv_20 + Jinv_01*Jinv_21 + Jinv_02*Jinv_22 - w[5][0]*(Jinv_00*Jinv_20 + Jinv_01*Jinv_21 + Jinv_02*Jinv_22));
    const double G13 = det*w[4][0]*w[5][0]*(Jinv_20*Jinv_20 + Jinv_21*Jinv_21 + Jinv_22*Jinv_22);
    const double G14 = det*w[4][0]*w[5][0]*(Jinv_00*Jinv_20 + Jinv_01*Jinv_21 + Jinv_02*Jinv_22);
    const double G15 =  - det;
    const double G16 =  - det*w[2][0]*(Jinv_20*Jinv_20 + Jinv_21*Jinv_21 + Jinv_22*Jinv_22);
    const double G17 = det*w[4][0]*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01 + Jinv_02*Jinv_02 - w[5][0]*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01 + Jinv_02*Jinv_02));
    const double G18 = det*w[4][0]*w[5][0]*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01 + Jinv_02*Jinv_02);
    const double G19 =  - 4*det*w[3][0];
    const double G20 = 6*det*w[3][0];
    const double G21 =  - 2*det*w[3][0];
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', True), ('ignore zero tables', True), ('non zero columns', True), ('remove zero terms', True), ('ignore ones', True)
    // Total number of operations to compute element tensor: 4671
    
    // Loop quadrature points for integral
    // Number of operations to compute element tensor for following IP loop = 4482
    for (unsigned int ip = 0; ip < 27; ip++)
    {
      
      // Function declarations
      double F0 = 0;
      double F1 = 0;
      double F2 = 0;
      double F3 = 0;
      double F4 = 0;
      double F5 = 0;
      double F6 = 0;
      double F7 = 0;
      double F8 = 0;
      double F9 = 0;
      double F10 = 0;
      double F11 = 0;
      
      // Total number of operations to compute function values = 36
      for (unsigned int r = 0; r < 2; r++)
      {
        F0 += FE0_C1_D001[ip][r]*w[0][nzc3[r]];
        F1 += FE0_C1_D001[ip][r]*w[0][nzc4[r]];
        F2 += FE0_C1_D001[ip][r]*w[0][nzc5[r]];
        F3 += FE0_C1_D001[ip][r]*w[1][nzc3[r]];
        F4 += FE0_C1_D001[ip][r]*w[1][nzc4[r]];
        F5 += FE0_C1_D001[ip][r]*w[1][nzc5[r]];
        F8 += FE0_C1_D001[ip][r]*w[0][nzc2[r]];
        F9 += FE0_C1_D001[ip][r]*w[0][nzc1[r]];
        F10 += FE0_C1_D001[ip][r]*w[0][nzc0[r]];
      }// end loop over 'r'
      
      // Total number of operations to compute function values = 24
      for (unsigned int r = 0; r < 4; r++)
      {
        F6 += FE0_C1[ip][r]*w[0][nzc6[r]];
        F7 += FE0_C1[ip][r]*w[1][nzc6[r]];
        F11 += FE0_C1[ip][r]*w[0][nzc7[r]];
      }// end loop over 'r'
      
      // Number of operations to compute ip constants: 66
      // Number of operations: 6
      const double Gip0 = W27[ip]*(F10*G0 + F8*G2 + F9*G1);
      
      // Number of operations: 12
      const double Gip1 = W27[ip]*(F0*G8 + F1*G5 + F2*G7 + F3*G6 + F4*G4 + F5*G3);
      
      // Number of operations: 6
      const double Gip2 = W27[ip]*(F10*G9 + F8*G10 + F9*G2);
      
      // Number of operations: 12
      const double Gip3 = W27[ip]*(F0*G14 + F1*G7 + F2*G13 + F3*G12 + F4*G3 + F5*G11);
      
      // Number of operations: 4
      const double Gip4 = W27[ip]*(F6*det + F7*G15);
      
      // Number of operations: 6
      const double Gip5 = W27[ip]*(F10*G16 + F8*G9 + F9*G0);
      
      // Number of operations: 12
      const double Gip6 = W27[ip]*(F0*G18 + F1*G8 + F2*G14 + F3*G17 + F4*G6 + F5*G12);
      
      // Number of operations: 8
      const double Gip7 = W27[ip]*(F11*det + F6*(G21 + F6*(G20 + F6*G19)));
      
      
      // Number of operations for primary indices = 24
      for (unsigned int j = 0; j < 2; j++)
      {
        // Number of operations to compute entry = 2
        A[nzc1[j]] += FE0_C1_D001[ip][j]*Gip0;
        // Number of operations to compute entry = 2
        A[nzc4[j]] += FE0_C1_D001[ip][j]*Gip1;
        // Number of operations to compute entry = 2
        A[nzc2[j]] += FE0_C1_D001[ip][j]*Gip2;
        // Number of operations to compute entry = 2
        A[nzc5[j]] += FE0_C1_D001[ip][j]*Gip3;
        // Number of operations to compute entry = 2
        A[nzc0[j]] += FE0_C1_D001[ip][j]*Gip5;
        // Number of operations to compute entry = 2
        A[nzc3[j]] += FE0_C1_D001[ip][j]*Gip6;
      }// end loop over 'j'
      
      // Number of operations for primary indices = 16
      for (unsigned int j = 0; j < 4; j++)
      {
        // Number of operations to compute entry = 2
        A[nzc7[j]] += FE0_C1[ip][j]*Gip4;
        // Number of operations to compute entry = 2
        A[nzc6[j]] += FE0_C1[ip][j]*Gip7;
      }// end loop over 'j'
    }// end loop over 'ip'
}

/// Constructor
cahnhilliard3d_1_cell_integral_0::cahnhilliard3d_1_cell_integral_0() : ufc::cell_integral()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_1_cell_integral_0::~cahnhilliard3d_1_cell_integral_0()
{
    // Do nothing
}

/// Tabulate the tensor for the contribution from a local cell
void cahnhilliard3d_1_cell_integral_0::tabulate_tensor(double* A,
                                    const double * const * w,
                                    const ufc::cell& c) const
{
    // Reset values of the element tensor block
    A[0] = 0;
    A[1] = 0;
    A[2] = 0;
    A[3] = 0;
    A[4] = 0;
    A[5] = 0;
    A[6] = 0;
    A[7] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c);
}

/// Constructor
cahnhilliard3d_form_1::cahnhilliard3d_form_1() : ufc::form()
{
    // Do nothing
}

/// Destructor
cahnhilliard3d_form_1::~cahnhilliard3d_form_1()
{
    // Do nothing
}

/// Return a string identifying the form
const char* cahnhilliard3d_form_1::signature() const
{
    return "Form([Integral(Sum(Sum(Product(Constant(Cell('tetrahedron', 1), 4), IndexSum(Product(Indexed(ComponentTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((Index(0),), {Index(0): 3})), MultiIndex((FixedIndex(0),), {})), MultiIndex((Index(0),), {Index(0): 3})), MultiIndex((Index(1),), {Index(1): 3})), Indexed(ComponentTensor(Sum(Product(Constant(Cell('tetrahedron', 1), 5), Indexed(SpatialDerivative(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((Index(2),), {Index(2): 3})), MultiIndex((FixedIndex(0),), {}))), Product(Indexed(SpatialDerivative(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((Index(2),), {Index(2): 3})), MultiIndex((FixedIndex(0),), {})), Sum(IntValue(1, (), (), {}), Product(IntValue(-1, (), (), {}), Constant(Cell('tetrahedron', 1), 5))))), MultiIndex((Index(2),), {Index(2): 3})), MultiIndex((Index(1),), {Index(1): 3}))), MultiIndex((Index(1),), {Index(1): 3}))), Sum(Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 2})), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2}))), Product(IntValue(-1, (), (), {}), Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 2})), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 1), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))))), Sum(Product(IntValue(-1, (), (), {}), Product(Constant(Cell('tetrahedron', 1), 2), IndexSum(Product(Indexed(ComponentTensor(Indexed(SpatialDerivative(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((Index(3),), {Index(3): 3})), MultiIndex((FixedIndex(1),), {})), MultiIndex((Index(3),), {Index(3): 3})), MultiIndex((Index(4),), {Index(4): 3})), Indexed(ComponentTensor(Indexed(SpatialDerivative(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((Index(5),), {Index(5): 3})), MultiIndex((FixedIndex(1),), {})), MultiIndex((Index(5),), {Index(5): 3})), MultiIndex((Index(4),), {Index(4): 3}))), MultiIndex((Index(4),), {Index(4): 3})))), Sum(Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(0),), {FixedIndex(0): 2}))), Product(IntValue(-1, (), (), {}), Product(Indexed(BasisFunction(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Product(Constant(Cell('tetrahedron', 1), 3), Sum(Product(IntValue(-1, (), (), {}), Product(Product(Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})), Product(IntValue(2, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))), Sum(IntValue(1, (), (), {}), Product(IntValue(-1, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))))), Product(Product(Product(IntValue(2, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2}))), Sum(IntValue(1, (), (), {}), Product(IntValue(-1, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2}))))), Sum(IntValue(1, (), (), {}), Product(IntValue(-1, (), (), {}), Indexed(Function(MixedElement(*[FiniteElement('Lagrange', Cell('tetrahedron', 1), 1), FiniteElement('Lagrange', Cell('tetrahedron', 1), 1)], **{'value_shape': (2,) }), 0), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2})))))))))))), Measure('cell', 0, None))])";
}

/// Return the rank of the global tensor (r)
unsigned int cahnhilliard3d_form_1::rank() const
{
    return 1;
}

/// Return the number of coefficients (n)
unsigned int cahnhilliard3d_form_1::num_coefficients() const
{
    return 6;
}

/// Return the number of cell integrals
unsigned int cahnhilliard3d_form_1::num_cell_integrals() const
{
    return 1;
}
  
/// Return the number of exterior facet integrals
unsigned int cahnhilliard3d_form_1::num_exterior_facet_integrals() const
{
    return 0;
}
  
/// Return the number of interior facet integrals
unsigned int cahnhilliard3d_form_1::num_interior_facet_integrals() const
{
    return 0;
}

/// Create a new finite element for argument function i
ufc::finite_element* cahnhilliard3d_form_1::create_finite_element(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_finite_element_0();
      break;
    case 1:
      return new cahnhilliard3d_1_finite_element_1();
      break;
    case 2:
      return new cahnhilliard3d_1_finite_element_2();
      break;
    case 3:
      return new cahnhilliard3d_1_finite_element_3();
      break;
    case 4:
      return new cahnhilliard3d_1_finite_element_4();
      break;
    case 5:
      return new cahnhilliard3d_1_finite_element_5();
      break;
    case 6:
      return new cahnhilliard3d_1_finite_element_6();
      break;
    }
    return 0;
}

/// Create a new dof map for argument function i
ufc::dof_map* cahnhilliard3d_form_1::create_dof_map(unsigned int i) const
{
    switch ( i )
    {
    case 0:
      return new cahnhilliard3d_1_dof_map_0();
      break;
    case 1:
      return new cahnhilliard3d_1_dof_map_1();
      break;
    case 2:
      return new cahnhilliard3d_1_dof_map_2();
      break;
    case 3:
      return new cahnhilliard3d_1_dof_map_3();
      break;
    case 4:
      return new cahnhilliard3d_1_dof_map_4();
      break;
    case 5:
      return new cahnhilliard3d_1_dof_map_5();
      break;
    case 6:
      return new cahnhilliard3d_1_dof_map_6();
      break;
    }
    return 0;
}

/// Create a new cell integral on sub domain i
ufc::cell_integral* cahnhilliard3d_form_1::create_cell_integral(unsigned int i) const
{
    return new cahnhilliard3d_1_cell_integral_0();
}

/// Create a new exterior facet integral on sub domain i
ufc::exterior_facet_integral* cahnhilliard3d_form_1::create_exterior_facet_integral(unsigned int i) const
{
    return 0;
}

/// Create a new interior facet integral on sub domain i
ufc::interior_facet_integral* cahnhilliard3d_form_1::create_interior_facet_integral(unsigned int i) const
{
    return 0;
}

