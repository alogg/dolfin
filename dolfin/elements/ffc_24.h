// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.5.1.

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
    values[0] = 0;
    values[1] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
    // Table(s) of coefficients
    const static double coefficients0[6][3] = \
    {{0.942809041582063, 0.577350269189626, -0.333333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.471404520791032, -0.577350269189626, -0.666666666666667},
    {0.471404520791032, 0.288675134594813, 0.833333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.942809041582064, 0.577350269189626, -0.333333333333333}};
    
    const static double coefficients1[6][3] = \
    {{-0.471404520791032, 0, -0.333333333333333},
    {0.942809041582063, 0, 0.666666666666667},
    {0.471404520791032, 0, 0.333333333333333},
    {-0.942809041582064, 0, -0.666666666666667},
    {-0.471404520791031, 0.866025403784439, 0.166666666666667},
    {-0.471404520791032, -0.866025403784439, 0.166666666666667}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    const double coeff1_0 = coefficients1[dof][0];
    const double coeff1_1 = coefficients1[dof][1];
    const double coeff1_2 = coefficients1[dof][2];
    
    // Compute value(s)
    const double tmp0_0 = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
    const double tmp0_1 = coeff1_0*basisvalue0 + coeff1_1*basisvalue1 + coeff1_2*basisvalue2;
    // Using contravariant Piola transform to map values back to the physical element
    values[0] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
    values[1] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
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
    for (unsigned int j = 0; j < 2*num_derivatives; j++)
      values[j] = 0;
    
    // Map degree of freedom to element degree of freedom
    const unsigned int dof = i;
    
    // Generate scalings
    const double scalings_y_0 = 1;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5*y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1;
    const double psitilde_a_1 = x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1;
    const double psitilde_bs_0_1 = 1.5*y + 0.5;
    const double psitilde_bs_1_0 = 1;
    
    // Compute basisvalues
    const double basisvalue0 = 0.707106781186548*psitilde_a_0*scalings_y_0*psitilde_bs_0_0;
    const double basisvalue1 = 1.73205080756888*psitilde_a_1*scalings_y_1*psitilde_bs_1_0;
    const double basisvalue2 = psitilde_a_0*scalings_y_0*psitilde_bs_0_1;
    
    // Table(s) of coefficients
    const static double coefficients0[6][3] = \
    {{0.942809041582063, 0.577350269189626, -0.333333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.471404520791032, -0.577350269189626, -0.666666666666667},
    {0.471404520791032, 0.288675134594813, 0.833333333333333},
    {-0.471404520791032, -0.288675134594813, 0.166666666666667},
    {0.942809041582064, 0.577350269189626, -0.333333333333333}};
    
    const static double coefficients1[6][3] = \
    {{-0.471404520791032, 0, -0.333333333333333},
    {0.942809041582063, 0, 0.666666666666667},
    {0.471404520791032, 0, 0.333333333333333},
    {-0.942809041582064, 0, -0.666666666666667},
    {-0.471404520791031, 0.866025403784439, 0.166666666666667},
    {-0.471404520791032, -0.866025403784439, 0.166666666666667}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    const static double dmats0[3][3] = \
    {{0, 0, 0},
    {4.89897948556636, 0, 0},
    {0, 0, 0}};
    
    const static double dmats1[3][3] = \
    {{0, 0, 0},
    {2.44948974278318, 0, 0},
    {4.24264068711928, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [2*num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    double coeff1_0 = 0;
    double coeff1_1 = 0;
    double coeff1_2 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    double new_coeff1_0 = 0;
    double new_coeff1_1 = 0;
    double new_coeff1_2 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
      new_coeff1_0 = coefficients1[dof][0];
      new_coeff1_1 = coefficients1[dof][1];
      new_coeff1_2 = coefficients1[dof][2];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
        coeff1_0 = new_coeff1_0;
        coeff1_1 = new_coeff1_1;
        coeff1_2 = new_coeff1_2;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2];
          new_coeff1_0 = coeff1_0*dmats0[0][0] + coeff1_1*dmats0[1][0] + coeff1_2*dmats0[2][0];
          new_coeff1_1 = coeff1_0*dmats0[0][1] + coeff1_1*dmats0[1][1] + coeff1_2*dmats0[2][1];
          new_coeff1_2 = coeff1_0*dmats0[0][2] + coeff1_1*dmats0[1][2] + coeff1_2*dmats0[2][2];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2];
          new_coeff1_0 = coeff1_0*dmats1[0][0] + coeff1_1*dmats1[1][0] + coeff1_2*dmats1[2][0];
          new_coeff1_1 = coeff1_0*dmats1[0][1] + coeff1_1*dmats1[1][1] + coeff1_2*dmats1[2][1];
          new_coeff1_2 = coeff1_0*dmats1[0][2] + coeff1_1*dmats1[1][2] + coeff1_2*dmats1[2][2];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      // Correct values by the contravariant Piola transform
      const double tmp0_0 = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2;
      const double tmp0_1 = new_coeff1_0*basisvalue0 + new_coeff1_1*basisvalue1 + new_coeff1_2*basisvalue2;
      derivatives[deriv_num] = (1.0/detJ)*(J_00*tmp0_0 + J_01*tmp0_1);
      derivatives[num_derivatives + deriv_num] = (1.0/detJ)*(J_10*tmp0_0 + J_11*tmp0_1);
    }
    
    // Transform derivatives back to physical element
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        values[row] += transform[row][col]*derivatives[col];
        values[num_derivatives + row] += transform[row][col]*derivatives[num_derivatives + col];
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
    vertex_values[2] = (1.0/detJ)*(dof_values[0]*2*J_00 + dof_values[1]*J_00 + dof_values[4]*(J_00 + J_01) + dof_values[5]*(2*J_00 - 2*J_01));
    vertex_values[4] = (1.0/detJ)*(dof_values[0]*J_01 + dof_values[1]*2*J_01 + dof_values[2]*(J_00 + J_01) + dof_values[3]*(2*J_00 - 2*J_01));
    vertex_values[1] = (1.0/detJ)*(dof_values[2]*2*J_10 + dof_values[3]*J_10 + dof_values[4]*(-2*J_11) + dof_values[5]*J_11);
    vertex_values[3] = (1.0/detJ)*(dof_values[0]*2*J_10 + dof_values[1]*J_10 + dof_values[4]*(J_10 + J_11) + dof_values[5]*(2*J_10 - 2*J_11));
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
    return 2;
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
