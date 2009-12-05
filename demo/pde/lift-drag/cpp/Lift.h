// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.7.1.
//
// Warning: This code was generated with the option '-l dolfin'
// and contains DOLFIN-specific wrappers that depend on DOLFIN.

#ifndef __LIFT_H
#define __LIFT_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>
    
/// This class defines the interface for a finite element.

class lift_0_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  lift_0_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~lift_0_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "FiniteElement('Lagrange', Cell('triangle', 1, Space(2)), 1)";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::triangle;
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
    return 3;
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
    static const double coefficients0[3][3] = \
    {{0.471404520791032, -0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0, 0.333333333333333}};
    
    // Extract relevant coefficients
    const double coeff0_0 = coefficients0[dof][0];
    const double coeff0_1 = coefficients0[dof][1];
    const double coeff0_2 = coefficients0[dof][2];
    
    // Compute value(s)
    *values = coeff0_0*basisvalue0 + coeff0_1*basisvalue1 + coeff0_2*basisvalue2;
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
    static const double coefficients0[3][3] = \
    {{0.471404520791032, -0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0.288675134594813, -0.166666666666667},
    {0.471404520791032, 0, 0.333333333333333}};
    
    // Interesting (new) part
    // Tables of derivatives of the polynomial base (transpose)
    static const double dmats0[3][3] = \
    {{0, 0, 0},
    {4.89897948556636, 0, 0},
    {0, 0, 0}};
    
    static const double dmats1[3][3] = \
    {{0, 0, 0},
    {2.44948974278318, 0, 0},
    {4.24264068711928, 0, 0}};
    
    // Compute reference derivatives
    // Declare pointer to array of derivatives on FIAT element
    double *derivatives = new double [num_derivatives];
    
    // Declare coefficients
    double coeff0_0 = 0;
    double coeff0_1 = 0;
    double coeff0_2 = 0;
    
    // Declare new coefficients
    double new_coeff0_0 = 0;
    double new_coeff0_1 = 0;
    double new_coeff0_2 = 0;
    
    // Loop possible derivatives
    for (unsigned int deriv_num = 0; deriv_num < num_derivatives; deriv_num++)
    {
      // Get values from coefficients array
      new_coeff0_0 = coefficients0[dof][0];
      new_coeff0_1 = coefficients0[dof][1];
      new_coeff0_2 = coefficients0[dof][2];
    
      // Loop derivative order
      for (unsigned int j = 0; j < n; j++)
      {
        // Update old coefficients
        coeff0_0 = new_coeff0_0;
        coeff0_1 = new_coeff0_1;
        coeff0_2 = new_coeff0_2;
    
        if(combinations[deriv_num][j] == 0)
        {
          new_coeff0_0 = coeff0_0*dmats0[0][0] + coeff0_1*dmats0[1][0] + coeff0_2*dmats0[2][0];
          new_coeff0_1 = coeff0_0*dmats0[0][1] + coeff0_1*dmats0[1][1] + coeff0_2*dmats0[2][1];
          new_coeff0_2 = coeff0_0*dmats0[0][2] + coeff0_1*dmats0[1][2] + coeff0_2*dmats0[2][2];
        }
        if(combinations[deriv_num][j] == 1)
        {
          new_coeff0_0 = coeff0_0*dmats1[0][0] + coeff0_1*dmats1[1][0] + coeff0_2*dmats1[2][0];
          new_coeff0_1 = coeff0_0*dmats1[0][1] + coeff0_1*dmats1[1][1] + coeff0_2*dmats1[2][1];
          new_coeff0_2 = coeff0_0*dmats1[0][2] + coeff0_1*dmats1[1][2] + coeff0_2*dmats1[2][2];
        }
    
      }
      // Compute derivatives on reference element as dot product of coefficients and basisvalues
      derivatives[deriv_num] = new_coeff0_0*basisvalue0 + new_coeff0_1*basisvalue1 + new_coeff0_2*basisvalue2;
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
    static const double X[3][1][2] = {{{0, 0}}, {{1, 0}}, {{0, 1}}};
    static const double W[3][1] = {{1}, {1}, {1}};
    static const double D[3][1][1] = {{{1}}, {{1}}, {{1}}};
    
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
    return new lift_0_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class lift_0_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  lift_0_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~lift_0_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for FiniteElement('Lagrange', Cell('triangle', 1, Space(2)), 1)";
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
      return false;
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
    __global_dimension = m.num_entities[0];
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
    return 3;
  }

  /// Return the maximum dimension of the local finite element function space
  virtual unsigned int max_local_dimension() const
  {
    return 3;
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
    dofs[0] = c.entity_indices[0][0];
    dofs[1] = c.entity_indices[0][1];
    dofs[2] = c.entity_indices[0][2];
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
      break;
    case 1:
      dofs[0] = 0;
      dofs[1] = 2;
      break;
    case 2:
      dofs[0] = 0;
      dofs[1] = 1;
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
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    return 1;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    return new lift_0_dof_map_0();
  }

};

/// This class defines the interface for the tabulation of the
/// exterior facet tensor corresponding to the local contribution to
/// a form from the integral over an exterior facet.

class lift_0_exterior_facet_integral_0_quadrature: public ufc::exterior_facet_integral
{
public:

  /// Constructor
  lift_0_exterior_facet_integral_0_quadrature() : ufc::exterior_facet_integral()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~lift_0_exterior_facet_integral_0_quadrature()
  {
    // Do nothing
  }

  /// Tabulate the tensor for the contribution from a local exterior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int facet) const
  {
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    
    // Compute determinant of Jacobian
    
    // Compute inverse of Jacobian
    
    // Vertices on edges
    static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    
    // Get vertices
    const unsigned int v0 = edge_vertices[facet][0];
    const unsigned int v1 = edge_vertices[facet][1];
    
    // Compute scale factor (length of edge scaled by length of reference interval)
    const double dx0 = x[v1][0] - x[v0][0];
    const double dx1 = x[v1][1] - x[v0][1];
    const double det = std::sqrt(dx0*dx0 + dx1*dx1);
    
    const bool direction = dx1*(x[facet][0] - x[v0][0]) - dx0*(x[facet][1] - x[v0][1]) < 0;
    
    // Compute facet normals from the facet scale factor constants
    const double n1 = direction ? -dx0 / det : dx0 / det;
    
    
    // Array of quadrature weights
    static const double W1 = 1;
    // Quadrature points on the UFC reference element: (0.5)
    
    // Value of basis functions at quadrature points.
    static const double FE0_f0[1][2] = \
    {{0.5, 0.5}};
    
    // Array of non-zero columns
    static const unsigned int nzc0[2] = {1, 2};
    // Array of non-zero columns
    static const unsigned int nzc1[2] = {0, 2};
    // Array of non-zero columns
    static const unsigned int nzc2[2] = {0, 1};
    
    // Number of operations to compute geometry constants: 2
    // Should be added to total operation count.
    const double G0 = W1*det*n1;
    
    // Compute element tensor using UFL quadrature representation
    // Optimisations: ('simplify expressions', True), ('ignore zero tables', True), ('non zero columns', True), ('remove zero terms', True), ('ignore ones', True)
    switch ( facet )
    {
    case 0:
      {
      // Total number of operations to compute element tensor (from this point): 6
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 6
      // Only 1 integration point, omitting IP loop.
      
      // Function declarations
      double F0 = 0;
      
      // Total number of operations to compute function values = 4
      for (unsigned int r = 0; r < 2; r++)
      {
        F0 += FE0_f0[0][r]*w[0][nzc0[r]];
      }// end loop over 'r'
      
      // Number of operations to compute ip constants: 1
      // Number of operations: 1
      const double Gip0 = F0*G0;
      
      
      // Number of operations for primary indices: 1
      // Number of operations to compute entry: 1
      A[0] += Gip0;
      }
      break;
    case 1:
      {
      // Total number of operations to compute element tensor (from this point): 6
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 6
      // Only 1 integration point, omitting IP loop.
      
      // Function declarations
      double F0 = 0;
      
      // Total number of operations to compute function values = 4
      for (unsigned int r = 0; r < 2; r++)
      {
        F0 += FE0_f0[0][r]*w[0][nzc1[r]];
      }// end loop over 'r'
      
      // Number of operations to compute ip constants: 1
      // Number of operations: 1
      const double Gip0 = F0*G0;
      
      
      // Number of operations for primary indices: 1
      // Number of operations to compute entry: 1
      A[0] += Gip0;
      }
      break;
    case 2:
      {
      // Total number of operations to compute element tensor (from this point): 6
      
      // Loop quadrature points for integral
      // Number of operations to compute element tensor for following IP loop = 6
      // Only 1 integration point, omitting IP loop.
      
      // Function declarations
      double F0 = 0;
      
      // Total number of operations to compute function values = 4
      for (unsigned int r = 0; r < 2; r++)
      {
        F0 += FE0_f0[0][r]*w[0][nzc2[r]];
      }// end loop over 'r'
      
      // Number of operations to compute ip constants: 1
      // Number of operations: 1
      const double Gip0 = F0*G0;
      
      
      // Number of operations for primary indices: 1
      // Number of operations to compute entry: 1
      A[0] += Gip0;
      }
      break;
    }
  }

};

/// This class defines the interface for the tabulation of the
/// exterior facet tensor corresponding to the local contribution to
/// a form from the integral over an exterior facet.

class lift_0_exterior_facet_integral_0: public ufc::exterior_facet_integral
{
private:

  lift_0_exterior_facet_integral_0_quadrature integral_0_quadrature;

public:

  /// Constructor
  lift_0_exterior_facet_integral_0() : ufc::exterior_facet_integral()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~lift_0_exterior_facet_integral_0()
  {
    // Do nothing
  }

  /// Tabulate the tensor for the contribution from a local exterior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int facet) const
  {
    // Reset values of the element tensor block
    A[0] = 0;
    
    // Add all contributions to element tensor
    integral_0_quadrature.tabulate_tensor(A, w, c, facet);
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

class lift_form_0: public ufc::form
{
public:

  /// Constructor
  lift_form_0() : ufc::form()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~lift_form_0()
  {
    // Do nothing
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
    return "Form([Integral(Product(Function(FiniteElement('Lagrange', Cell('triangle', 1, Space(2)), 1), 0), Indexed(FacetNormal(Cell('triangle', 1, Space(2))), MultiIndex((FixedIndex(1),), {FixedIndex(1): 2}))), Measure('exterior_facet', 0, None))])";
  }

  /// Return the rank of the global tensor (r)
  virtual unsigned int rank() const
  {
    return 0;
  }

  /// Return the number of coefficients (n)
  virtual unsigned int num_coefficients() const
  {
    return 1;
  }

  /// Return the number of cell integrals
  virtual unsigned int num_cell_integrals() const
  {
    return 0;
  }

  /// Return the number of exterior facet integrals
  virtual unsigned int num_exterior_facet_integrals() const
  {
    return 1;
  }

  /// Return the number of interior facet integrals
  virtual unsigned int num_interior_facet_integrals() const
  {
    return 0;
  }

  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(unsigned int i) const
  {
    return new lift_0_finite_element_0();
  }

  /// Create a new dof map for argument function i
  virtual ufc::dof_map* create_dof_map(unsigned int i) const
  {
    return new lift_0_dof_map_0();
  }

  /// Create a new cell integral on sub domain i
  virtual ufc::cell_integral* create_cell_integral(unsigned int i) const
  {
    return 0;
  }

  /// Create a new exterior facet integral on sub domain i
  virtual ufc::exterior_facet_integral* create_exterior_facet_integral(unsigned int i) const
  {
    return new lift_0_exterior_facet_integral_0();
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

namespace Lift
{

class CoefficientSpace_p: public dolfin::FunctionSpace
{
public:

  CoefficientSpace_p(const dolfin::Mesh& mesh):
      dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new lift_0_finite_element_0()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new lift_0_dof_map_0()), mesh)))
  {
    // Do nothing
  }

  CoefficientSpace_p(dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new lift_0_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new lift_0_dof_map_0()), mesh)))
  {
    // Do nothing
  }

  CoefficientSpace_p(boost::shared_ptr<dolfin::Mesh> mesh):
      dolfin::FunctionSpace(mesh,
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new lift_0_finite_element_0()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new lift_0_dof_map_0()), *mesh)))
  {
      // Do nothing
  }

  CoefficientSpace_p(boost::shared_ptr<const dolfin::Mesh> mesh):
      dolfin::FunctionSpace(mesh,
                            boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new lift_0_finite_element_0()))),
                            boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dof_map>(new lift_0_dof_map_0()), *mesh)))
  {
      // Do nothing
  }


  ~CoefficientSpace_p()
  {
  }

};

typedef CoefficientSpace_p Form_0_FunctionSpace_0;

class Form_0: public dolfin::Form
{
public:

  // Constructor
  Form_0(const dolfin::Mesh& mesh):
    dolfin::Form(0, 1), p(*this, 0)
  {
    _mesh = reference_to_no_delete_pointer(mesh);
    _ufc_form = boost::shared_ptr<const ufc::form>(new lift_form_0());
  }

  // Constructor
  Form_0(const dolfin::Mesh& mesh, const dolfin::GenericFunction& p):
    dolfin::Form(0, 1), p(*this, 0)
  {
    _mesh = reference_to_no_delete_pointer(mesh);
    this->p = p;

    _ufc_form = boost::shared_ptr<const ufc::form>(new lift_form_0());
  }

  // Constructor
  Form_0(const dolfin::Mesh& mesh, boost::shared_ptr<const dolfin::GenericFunction> p):
    dolfin::Form(0, 1), p(*this, 0)
  {
    _mesh = reference_to_no_delete_pointer(mesh);
    this->p = *p;

    _ufc_form = boost::shared_ptr<const ufc::form>(new lift_form_0());
  }

  // Constructor
  Form_0(boost::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::Form(0, 1), p(*this, 0)
  {
    _mesh = mesh;
    _ufc_form = boost::shared_ptr<const ufc::form>(new lift_form_0());
  }

  // Constructor
  Form_0(boost::shared_ptr<const dolfin::Mesh> mesh, const dolfin::GenericFunction& p):
    dolfin::Form(0, 1), p(*this, 0)
  {
    _mesh = mesh;
    this->p = p;

    _ufc_form = boost::shared_ptr<const ufc::form>(new lift_form_0());
  }

  // Constructor
  Form_0(boost::shared_ptr<const dolfin::Mesh> mesh, boost::shared_ptr<const dolfin::GenericFunction> p):
    dolfin::Form(0, 1), p(*this, 0)
  {
    _mesh = mesh;
    this->p = *p;

    _ufc_form = boost::shared_ptr<const ufc::form>(new lift_form_0());
  }

  // Destructor
  ~Form_0()
  {}

  /// Return the number of the coefficient with this name
  virtual dolfin::uint coefficient_number(const std::string& name) const
  {
    if (name == "p")
      return 0;

    dolfin::error("Invalid coefficient.");
    return 0;
  }

  /// Return the name of the coefficient with this number
  virtual std::string coefficient_name(dolfin::uint i) const
  {
    switch (i)
    {
    case 0:
      return "p";
    }

    dolfin::error("Invalid coefficient.");
    return "unnamed";
  }

  // Typedefs
  typedef Form_0_FunctionSpace_0 CoefficientSpace_p;

  // Coefficients
  dolfin::CoefficientAssigner p;
};

// Class typedefs
typedef Form_0 Functional;

} // namespace Lift

#endif
