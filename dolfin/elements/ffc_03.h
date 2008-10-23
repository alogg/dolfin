// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.5.1.

#ifndef __FFC_03_H
#define __FFC_03_H

#include <cmath>
#include <stdexcept>
#include <ufc.h>

/// This class defines the interface for a finite element.

class ffc_03_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  ffc_03_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~ffc_03_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Lagrange finite element of degree 2 on a triangle";
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
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    const static double coefficients0[6][6] = \
    {{0, -0.173205080756888, -0.1, 0.121716123890037, 0.0942809041582064, 0.0544331053951817},
    {0, 0.173205080756888, -0.1, 0.121716123890037, -0.0942809041582063, 0.0544331053951818},
    {0, 0, 0.2, 0, 0, 0.163299316185545},
    {0.471404520791032, 0.23094010767585, 0.133333333333333, 0, 0.188561808316413, -0.163299316185545},
    {0.471404520791032, -0.23094010767585, 0.133333333333333, 0, -0.188561808316413, -0.163299316185545},
    {0.471404520791032, 0, -0.266666666666667, -0.243432247780074, 0, 0.0544331053951817}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff0_3 = coefficients0[dof][3];
    const double coeff0_4 = coefficients0[dof][4];
    const double coeff0_5 = coefficients0[dof][5];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2 + coeff0_3*basisvalue3 + coeff0_4*basisvalue4 + coeff0_5*basisvalue5;
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
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    const double basisvalue3 = 2.73861278752583*psitilde_a_2*scalings_y_2*psitilde_bs_2_0;
    const double basisvalue4 = 2.12132034355964*psitilde_a_1*scalings_y_1*psitilde_bs_1_1;
    const double basisvalue5 = 1.22474487139159*psitilde_a_0*scalings_y_0*psitilde_bs_0_2;
    
    // Table(s) of coefficients
    const static double coefficients0[6][6] = \
    {{0, -0.173205080756888, -0.1, 0.121716123890037, 0.0942809041582064, 0.0544331053951817},
    {0, 0.173205080756888, -0.1, 0.121716123890037, -0.0942809041582063, 0.0544331053951818},
    {0, 0, 0.2, 0, 0, 0.163299316185545},
    {0.471404520791032, 0.23094010767585, 0.133333333333333, 0, 0.188561808316413, -0.163299316185545},
    {0.471404520791032, -0.23094010767585, 0.133333333333333, 0, -0.188561808316413, -0.163299316185545},
    {0.471404520791032, 0, -0.266666666666667, -0.243432247780074, 0, 0.0544331053951817}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {4.89897948556635, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 9.48683298050514, 0, 0, 0, 0},
    {4, 0, 7.07106781186548, 0, 0, 0},
    {0, 0, 0, 0, 0, 0}};
    
    const static double dmats1[6][6] = \
    {{0, 0, 0, 0, 0, 0},
    {2.44948974278318, 0, 0, 0, 0, 0},
    {4.24264068711928, 0, 0, 0, 0, 0},
    {2.58198889747161, 4.74341649025257, -0.912870929175277, 0, 0, 0},
    {2, 6.12372435695795, 3.53553390593274, 0, 0, 0},
    {-2.3094010767585, 0, 8.16496580927726, 0, 0, 0}};
    
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
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff0_3 = 0;
    double new_coeff0_4 = 0;
    double new_coeff0_5 = 0;
    
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
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0] + coeff0_3*dmats0[3][0] + coeff0_4*dmats0[4][0] + coeff0_5*dmats0[5][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1] + coeff0_3*dmats0[3][1] + coeff0_4*dmats0[4][1] + coeff0_5*dmats0[5][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2] + coeff0_3*dmats0[3][2] + coeff0_4*dmats0[4][2] + coeff0_5*dmats0[5][2];
          new_coeff0_3 = coeff0_0*dmats0[0][3] + coeff0_1*dmats0[1][3] + coeff0_2*dmats0[2][3] + coeff0_3*dmats0[3][3] + coeff0_4*dmats0[4][3] + coeff0_5*dmats0[5][3];
          new_coeff0_4 = coeff0_0*dmats0[0][4] + coeff0_1*dmats0[1][4] + coeff0_2*dmats0[2][4] + coeff0_3*dmats0[3][4] + coeff0_4*dmats0[4][4] + coeff0_5*dmats0[5][4];
          new_coeff0_5 = coeff0_0*dmats0[0][5] + coeff0_1*dmats0[1][5] + coeff0_2*dmats0[2][5] + coeff0_3*dmats0[3][5] + coeff0_4*dmats0[4][5] + coeff0_5*dmats0[5][5];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0] + coeff0_3*dmats1[3][0] + coeff0_4*dmats1[4][0] + coeff0_5*dmats1[5][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1] + coeff0_3*dmats1[3][1] + coeff0_4*dmats1[4][1] + coeff0_5*dmats1[5][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2] + coeff0_3*dmats1[3][2] + coeff0_4*dmats1[4][2] + coeff0_5*dmats1[5][2];
          new_coeff0_3 = coeff0_0*dmats1[0][3] + coeff0_1*dmats1[1][3] + coeff0_2*dmats1[2][3] + coeff0_3*dmats1[3][3] + coeff0_4*dmats1[4][3] + coeff0_5*dmats1[5][3];
          new_coeff0_4 = coeff0_0*dmats1[0][4] + coeff0_1*dmats1[1][4] + coeff0_2*dmats1[2][4] + coeff0_3*dmats1[3][4] + coeff0_4*dmats1[4][4] + coeff0_5*dmats1[5][4];
          new_coeff0_5 = coeff0_0*dmats1[0][5] + coeff0_1*dmats1[1][5] + coeff0_2*dmats1[2][5] + coeff0_3*dmats1[3][5] + coeff0_4*dmats1[4][5] + coeff0_5*dmats1[5][5];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2 + new_coeff0_3*basisvalue3 + new_coeff0_4*basisvalue4 + new_coeff0_5*basisvalue5;
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
    const static double X[6][1][2] = {{{0, 0}}, {{1, 0}}, {{0, 1}}, {{0.5, 0.5}}, {{0, 0.5}}, {{0.5, 0}}};
    const static double W[6][1] = {{1}, {1}, {1}, {1}, {1}, {1}};
    const static double D[6][1][1] = {{{1}}, {{1}}, {{1}}, {{1}}, {{1}}, {{1}}};
    
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
    return new ffc_03_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class ffc_03_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  ffc_03_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~ffc_03_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Lagrange finite element of degree 2 on a triangle";
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
    return 6;
  }

  // Return the geometric dimension of the coordinates this dof map provides
  virtual unsigned int geometric_dimension() const
  {
    return 2;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 3;
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
    dofs[3] = offset + c.entity_indices[1][0];
    dofs[4] = offset + c.entity_indices[1][1];
    dofs[5] = offset + c.entity_indices[1][2];
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
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      dofs[2] = 4;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
      dofs[2] = 5;
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
    coordinates[3][0] = 0.5*x[1][0] + 0.5*x[2][0];
    coordinates[3][1] = 0.5*x[1][1] + 0.5*x[2][1];
    coordinates[4][0] = 0.5*x[0][0] + 0.5*x[2][0];
    coordinates[4][1] = 0.5*x[0][1] + 0.5*x[2][1];
    coordinates[5][0] = 0.5*x[0][0] + 0.5*x[1][0];
    coordinates[5][1] = 0.5*x[0][1] + 0.5*x[1][1];
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new ffc_03_dof_map_0();
  }

};

#endif
