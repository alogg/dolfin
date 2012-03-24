// This code conforms with the UFC specification version 2.0.5
// and was automatically generated by FFC version 1.0.0+.
//
// This code was generated with the option '-l dolfin' and
// contains DOLFIN-specific wrappers that depend on DOLFIN.
// 
// This code was generated with the following parameters:
// 
//   cache_dir:                      ''
//   convert_exceptions_to_warnings: False
//   cpp_optimize:                   False
//   cpp_optimize_flags:             '-O2'
//   epsilon:                        1e-14
//   error_control:                  False
//   form_postfix:                   True
//   format:                         'dolfin'
//   log_level:                      10
//   log_prefix:                     ''
//   no_ferari:                      True
//   optimize:                       True
//   output_dir:                     '.'
//   precision:                      15
//   quadrature_degree:              'auto'
//   quadrature_rule:                'auto'
//   representation:                 'auto'
//   split:                          False
//   swig_binary:                    'swig'
//   swig_path:                      ''

#ifndef __POISSON1D_H
#define __POISSON1D_H

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <ufc.h>

/// This class defines the interface for a finite element.

class poisson1d_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  poisson1d_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~poisson1d_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)";
  }

  /// Return the cell shape
  virtual ufc::shape cell_shape() const
  {
    return ufc::interval;
  }

  /// Return the topological dimension of the cell shape
  virtual unsigned int topological_dimension() const
  {
    return 1;
  }

  /// Return the geometric dimension of the cell shape
  virtual unsigned int geometric_dimension() const
  {
    return 1;
  }

  /// Return the dimension of the finite element function space
  virtual unsigned int space_dimension() const
  {
    return 2;
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
    
    // Compute determinant of Jacobian
    
    // Compute inverse of Jacobian
    
    // Get coordinates and map to the reference (FIAT) element
    double X = (2.0*coordinates[0] - x[0][0] - x[1][0]) / J_00;
    
    // Reset values.
    *values = 0.0;
    switch (i)
    {
    case 0:
      {
        
      // Array of basisvalues.
      double basisvalues[2] = {0.0, 0.0};
      
      // Declare helper variables.
      
      // Compute basisvalues.
      basisvalues[0] = 1.0;
      basisvalues[1] = X;
      for (unsigned int r = 0; r < 2; r++)
      {
        basisvalues[r] *= std::sqrt((0.5 + r));
      }// end loop over 'r'
      
      // Table(s) of coefficients.
      static const double coefficients0[2] = \
      {0.707106781186547, -0.408248290463863};
      
      // Compute value(s).
      for (unsigned int r = 0; r < 2; r++)
      {
        *values += coefficients0[r]*basisvalues[r];
      }// end loop over 'r'
        break;
      }
    case 1:
      {
        
      // Array of basisvalues.
      double basisvalues[2] = {0.0, 0.0};
      
      // Declare helper variables.
      
      // Compute basisvalues.
      basisvalues[0] = 1.0;
      basisvalues[1] = X;
      for (unsigned int r = 0; r < 2; r++)
      {
        basisvalues[r] *= std::sqrt((0.5 + r));
      }// end loop over 'r'
      
      // Table(s) of coefficients.
      static const double coefficients0[2] = \
      {0.707106781186547, 0.408248290463863};
      
      // Compute value(s).
      for (unsigned int r = 0; r < 2; r++)
      {
        *values += coefficients0[r]*basisvalues[r];
      }// end loop over 'r'
        break;
      }
    }
    
  }

  /// Evaluate all basis functions at given point in cell
  virtual void evaluate_basis_all(double* values,
                                  const double* coordinates,
                                  const ufc::cell& c) const
  {
    // Helper variable to hold values of a single dof.
    double dof_values = 0.0;
    
    // Loop dofs and call evaluate_basis.
    for (unsigned int r = 0; r < 2; r++)
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
    
    // Compute determinant of Jacobian
    const double detJ =  J_00;
    
    // Compute inverse of Jacobian
    const double K_00 =  1.0 / detJ;
    
    // Get coordinates and map to the reference (FIAT) element
    double X = (2.0*coordinates[0] - x[0][0] - x[1][0]) / J_00;
    
    // Compute number of derivatives.
    unsigned int num_derivatives = 1;
    for (unsigned int r = 0; r < n; r++)
    {
      num_derivatives *= 1;
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
          if (combinations[row][col] + 1 > 0)
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
    const double Jinv[1][1] =  {{K_00}};
    
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
      values[r] = 0.0;
    }// end loop over 'r'
    
    switch (i)
    {
    case 0:
      {
        
      // Array of basisvalues.
      double basisvalues[2] = {0.0, 0.0};
      
      // Declare helper variables.
      
      // Compute basisvalues.
      basisvalues[0] = 1.0;
      basisvalues[1] = X;
      for (unsigned int r = 0; r < 2; r++)
      {
        basisvalues[r] *= std::sqrt((0.5 + r));
      }// end loop over 'r'
      
      // Table(s) of coefficients.
      static const double coefficients0[2] = \
      {0.707106781186547, -0.408248290463863};
      
      // Tables of derivatives of the polynomial base (transpose).
      static const double dmats0[2][2] = \
      {{0.0, 0.0},
      {3.46410161513775, 0.0}};
      
      // Compute reference derivatives.
      // Declare pointer to array of derivatives on FIAT element.
      double *derivatives = new double[num_derivatives];
      for (unsigned int r = 0; r < num_derivatives; r++)
      {
        derivatives[r] = 0.0;
      }// end loop over 'r'
      
      // Declare derivative matrix (of polynomial basis).
      double dmats[2][2] = \
      {{1.0, 0.0},
      {0.0, 1.0}};
      
      // Declare (auxiliary) derivative matrix (of polynomial basis).
      double dmats_old[2][2] = \
      {{1.0, 0.0},
      {0.0, 1.0}};
      
      // Loop possible derivatives.
      for (unsigned int r = 0; r < num_derivatives; r++)
      {
        // Resetting dmats values to compute next derivative.
        for (unsigned int t = 0; t < 2; t++)
        {
          for (unsigned int u = 0; u < 2; u++)
          {
            dmats[t][u] = 0.0;
            if (t == u)
            {
            dmats[t][u] = 1.0;
            }
            
          }// end loop over 'u'
        }// end loop over 't'
        
        // Looping derivative order to generate dmats.
        for (unsigned int s = 0; s < n; s++)
        {
          // Updating dmats_old with new values and resetting dmats.
          for (unsigned int t = 0; t < 2; t++)
          {
            for (unsigned int u = 0; u < 2; u++)
            {
              dmats_old[t][u] = dmats[t][u];
              dmats[t][u] = 0.0;
            }// end loop over 'u'
          }// end loop over 't'
          
          // Update dmats using an inner product.
          if (combinations[r][s] == 0)
          {
          for (unsigned int t = 0; t < 2; t++)
          {
            for (unsigned int u = 0; u < 2; u++)
            {
              for (unsigned int tu = 0; tu < 2; tu++)
              {
                dmats[t][u] += dmats0[t][tu]*dmats_old[tu][u];
              }// end loop over 'tu'
            }// end loop over 'u'
          }// end loop over 't'
          }
          
        }// end loop over 's'
        for (unsigned int s = 0; s < 2; s++)
        {
          for (unsigned int t = 0; t < 2; t++)
          {
            derivatives[r] += coefficients0[s]*dmats[s][t]*basisvalues[t];
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
        break;
      }
    case 1:
      {
        
      // Array of basisvalues.
      double basisvalues[2] = {0.0, 0.0};
      
      // Declare helper variables.
      
      // Compute basisvalues.
      basisvalues[0] = 1.0;
      basisvalues[1] = X;
      for (unsigned int r = 0; r < 2; r++)
      {
        basisvalues[r] *= std::sqrt((0.5 + r));
      }// end loop over 'r'
      
      // Table(s) of coefficients.
      static const double coefficients0[2] = \
      {0.707106781186547, 0.408248290463863};
      
      // Tables of derivatives of the polynomial base (transpose).
      static const double dmats0[2][2] = \
      {{0.0, 0.0},
      {3.46410161513775, 0.0}};
      
      // Compute reference derivatives.
      // Declare pointer to array of derivatives on FIAT element.
      double *derivatives = new double[num_derivatives];
      for (unsigned int r = 0; r < num_derivatives; r++)
      {
        derivatives[r] = 0.0;
      }// end loop over 'r'
      
      // Declare derivative matrix (of polynomial basis).
      double dmats[2][2] = \
      {{1.0, 0.0},
      {0.0, 1.0}};
      
      // Declare (auxiliary) derivative matrix (of polynomial basis).
      double dmats_old[2][2] = \
      {{1.0, 0.0},
      {0.0, 1.0}};
      
      // Loop possible derivatives.
      for (unsigned int r = 0; r < num_derivatives; r++)
      {
        // Resetting dmats values to compute next derivative.
        for (unsigned int t = 0; t < 2; t++)
        {
          for (unsigned int u = 0; u < 2; u++)
          {
            dmats[t][u] = 0.0;
            if (t == u)
            {
            dmats[t][u] = 1.0;
            }
            
          }// end loop over 'u'
        }// end loop over 't'
        
        // Looping derivative order to generate dmats.
        for (unsigned int s = 0; s < n; s++)
        {
          // Updating dmats_old with new values and resetting dmats.
          for (unsigned int t = 0; t < 2; t++)
          {
            for (unsigned int u = 0; u < 2; u++)
            {
              dmats_old[t][u] = dmats[t][u];
              dmats[t][u] = 0.0;
            }// end loop over 'u'
          }// end loop over 't'
          
          // Update dmats using an inner product.
          if (combinations[r][s] == 0)
          {
          for (unsigned int t = 0; t < 2; t++)
          {
            for (unsigned int u = 0; u < 2; u++)
            {
              for (unsigned int tu = 0; tu < 2; tu++)
              {
                dmats[t][u] += dmats0[t][tu]*dmats_old[tu][u];
              }// end loop over 'tu'
            }// end loop over 'u'
          }// end loop over 't'
          }
          
        }// end loop over 's'
        for (unsigned int s = 0; s < 2; s++)
        {
          for (unsigned int t = 0; t < 2; t++)
          {
            derivatives[r] += coefficients0[s]*dmats[s][t]*basisvalues[t];
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
        break;
      }
    }
    
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
      num_derivatives *= 1;
    }// end loop over 'r'
    
    // Helper variable to hold values of a single dof.
    double *dof_values = new double[num_derivatives];
    for (unsigned int r = 0; r < num_derivatives; r++)
    {
      dof_values[r] = 0.0;
    }// end loop over 'r'
    
    // Loop dofs and call evaluate_basis_derivatives.
    for (unsigned int r = 0; r < 2; r++)
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
    double y[1];
    const double * const * x = c.coordinates;
    switch (i)
    {
    case 0:
      {
        y[0] = x[0][0];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    case 1:
      {
        y[0] = x[1][0];
      f.evaluate(vals, y, c);
      return vals[0];
        break;
      }
    }
    
    return 0.0;
  }

  /// Evaluate linear functionals for all dofs on the function f
  virtual void evaluate_dofs(double* values,
                             const ufc::function& f,
                             const ufc::cell& c) const
  {
    // Declare variables for result of evaluation.
    double vals[1];
    
    // Declare variable for physical coordinates.
    double y[1];
    const double * const * x = c.coordinates;
    y[0] = x[0][0];
    f.evaluate(vals, y, c);
    values[0] = vals[0];
    y[0] = x[1][0];
    f.evaluate(vals, y, c);
    values[1] = vals[0];
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
    // Evaluate function and change variables
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
  }

  /// Map coordinate xhat from reference cell to coordinate x in cell
  virtual void map_from_reference_cell(double* x,
                                       const double* xhat,
                                       const ufc::cell& c) const
  {
    throw std::runtime_error("map_from_reference_cell not yet implemented (introduced in UFC 2.0).");
  }

  /// Map from coordinate x in cell to coordinate xhat in reference cell
  virtual void map_to_reference_cell(double* xhat,
                                     const double* x,
                                     const ufc::cell& c) const
  {
    throw std::runtime_error("map_to_reference_cell not yet implemented (introduced in UFC 2.0).");
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

  /// Create a new class instance
  virtual ufc::finite_element* create() const
  {
    return new poisson1d_finite_element_0();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class poisson1d_dofmap_0: public ufc::dofmap
{
private:

  unsigned int _global_dimension;
public:

  /// Constructor
  poisson1d_dofmap_0() : ufc::dofmap()
  {
    _global_dimension = 0;
  }

  /// Destructor
  virtual ~poisson1d_dofmap_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dofmap
  virtual const char* signature() const
  {
    return "FFC dofmap for FiniteElement('Lagrange', Cell('interval', Space(1)), 1, None)";
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
    }
    
    return false;
  }

  /// Initialize dofmap for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    _global_dimension = m.num_entities[0];
    return false;
  }

  /// Initialize dofmap for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c)
  {
    // Do nothing
  }

  /// Finish initialization of dofmap for cells
  virtual void init_cell_finalize()
  {
    // Do nothing
  }

  /// Return the topological dimension of the associated cell shape
  virtual unsigned int topological_dimension() const
  {
    return 1;
  }

  /// Return the geometric dimension of the associated cell shape
  virtual unsigned int geometric_dimension() const
  {
    return 1;
  }

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const
  {
    return _global_dimension;
  }

  /// Return the dimension of the local finite element function space for a cell
  virtual unsigned int local_dimension(const ufc::cell& c) const
  {
    return 2;
  }

  /// Return the maximum dimension of the local finite element function space
  virtual unsigned int max_local_dimension() const
  {
    return 2;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 1;
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
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   unsigned int facet) const
  {
    switch (facet)
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
    }
    
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(unsigned int* dofs,
                                    unsigned int d, unsigned int i) const
  {
    if (d > 1)
    {
    throw std::runtime_error("d is larger than dimension (1)");
    }
    
    switch (d)
    {
    case 0:
      {
        if (i > 1)
      {
      throw std::runtime_error("i is larger than number of entities (1)");
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
      }
      
        break;
      }
    case 1:
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
    coordinates[1][0] = x[1][0];
  }

  /// Return the number of sub dofmaps (for a mixed element)
  virtual unsigned int num_sub_dofmaps() const
  {
    return 0;
  }

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  virtual ufc::dofmap* create_sub_dofmap(unsigned int i) const
  {
    return 0;
  }

  /// Create a new class instance
  virtual ufc::dofmap* create() const
  {
    return new poisson1d_dofmap_0();
  }

};

/// This class defines the interface for the tabulation of the cell
/// tensor corresponding to the local contribution to a form from
/// the integral over a cell.

class poisson1d_cell_integral_0_0: public ufc::cell_integral
{
public:

  /// Constructor
  poisson1d_cell_integral_0_0() : ufc::cell_integral()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~poisson1d_cell_integral_0_0()
  {
    // Do nothing
  }

  /// Tabulate the tensor for the contribution from a local cell
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c) const
  {
    // Number of operations (multiply-add pairs) for Jacobian data:      7
    // Number of operations (multiply-add pairs) for geometry tensor:    1
    // Number of operations (multiply-add pairs) for tensor contraction: 0
    // Total number of operations (multiply-add pairs):                  8
    
    // Extract vertex coordinates
    const double * const * x = c.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J_00 = x[1][0] - x[0][0];
    
    // Compute determinant of Jacobian
    const double detJ =  J_00;
    
    // Compute inverse of Jacobian
    const double K_00 =  1.0 / detJ;
    
    // Set scale factor
    const double det = std::abs(detJ);
    
    // Compute geometry tensor
    const double G0_0_0 = det*K_00*K_00*(1.0);
    
    // Compute element tensor
    A[0] = G0_0_0;
    A[1] = -G0_0_0;
    A[2] = -G0_0_0;
    A[3] = G0_0_0;
  }

  /// Tabulate the tensor for the contribution from a local cell
  /// using the specified reference cell quadrature points/weights
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c,
                               unsigned int num_quadrature_points,
                               const double * const * quadrature_points,
                               const double* quadrature_weights) const
  {
    throw std::runtime_error("Quadrature version of tabulate_tensor not available when using the FFC tensor representation.");
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

class poisson1d_form_0: public ufc::form
{
public:

  /// Constructor
  poisson1d_form_0() : ufc::form()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~poisson1d_form_0()
  {
    // Do nothing
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
    return "e3c5e1f4caf29c3fcdaadf691af2167d1a0cec3f6dc24782fb213e4f7a6a2730e5d1131946f07146e1d7a24588879e9a691f1144778c36b56b955f17472717bc";
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

  /// Return the number of cell domains
  virtual unsigned int num_cell_domains() const
  {
    return 1;
  }

  /// Return the number of exterior facet domains
  virtual unsigned int num_exterior_facet_domains() const
  {
    return 0;
  }

  /// Return the number of interior facet domains
  virtual unsigned int num_interior_facet_domains() const
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
        return new poisson1d_finite_element_0();
        break;
      }
    case 1:
      {
        return new poisson1d_finite_element_0();
        break;
      }
    }
    
    return 0;
  }

  /// Create a new dofmap for argument function i
  virtual ufc::dofmap* create_dofmap(unsigned int i) const
  {
    switch (i)
    {
    case 0:
      {
        return new poisson1d_dofmap_0();
        break;
      }
    case 1:
      {
        return new poisson1d_dofmap_0();
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
        return new poisson1d_cell_integral_0_0();
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
#include <dolfin/adaptivity/ErrorControl.h>
#include <dolfin/adaptivity/GoalFunctional.h>

namespace Poisson1D
{

class Form_0_FunctionSpace_0: public dolfin::FunctionSpace
{
public:

  Form_0_FunctionSpace_0(const dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_0(dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_0(boost::shared_ptr<dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), *mesh)))
  {
      // Do nothing
  }

  Form_0_FunctionSpace_0(boost::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), *mesh)))
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
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_1(dolfin::Mesh& mesh):
    dolfin::FunctionSpace(dolfin::reference_to_no_delete_pointer(mesh),
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), mesh)))
  {
    // Do nothing
  }

  Form_0_FunctionSpace_1(boost::shared_ptr<dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), *mesh)))
  {
      // Do nothing
  }

  Form_0_FunctionSpace_1(boost::shared_ptr<const dolfin::Mesh> mesh):
    dolfin::FunctionSpace(mesh,
                          boost::shared_ptr<const dolfin::FiniteElement>(new dolfin::FiniteElement(boost::shared_ptr<ufc::finite_element>(new poisson1d_finite_element_0()))),
                          boost::shared_ptr<const dolfin::DofMap>(new dolfin::DofMap(boost::shared_ptr<ufc::dofmap>(new poisson1d_dofmap_0()), *mesh)))
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
  Form_0(const dolfin::FunctionSpace& V1, const dolfin::FunctionSpace& V0):
    dolfin::Form(2, 0)
  {
    _function_spaces[0] = reference_to_no_delete_pointer(V0);
    _function_spaces[1] = reference_to_no_delete_pointer(V1);

    _ufc_form = boost::shared_ptr<const ufc::form>(new poisson1d_form_0());
  }

  // Constructor
  Form_0(boost::shared_ptr<const dolfin::FunctionSpace> V1, boost::shared_ptr<const dolfin::FunctionSpace> V0):
    dolfin::Form(2, 0)
  {
    _function_spaces[0] = V0;
    _function_spaces[1] = V1;

    _ufc_form = boost::shared_ptr<const ufc::form>(new poisson1d_form_0());
  }

  // Destructor
  ~Form_0()
  {}

  /// Return the number of the coefficient with this name
  virtual dolfin::uint coefficient_number(const std::string& name) const
  {

    dolfin::dolfin_error("generated code for class Form",
                         "access coefficient data",
                         "There are no coefficients");
    return 0;
  }

  /// Return the name of the coefficient with this number
  virtual std::string coefficient_name(dolfin::uint i) const
  {

    dolfin::dolfin_error("generated code for class Form",
                         "access coefficient data",
                         "There are no coefficients");
    return "unnamed";
  }

  // Typedefs
  typedef Form_0_FunctionSpace_0 TestSpace;
  typedef Form_0_FunctionSpace_1 TrialSpace;

  // Coefficients
};

// Class typedefs
typedef Form_0 BilinearForm;
typedef Form_0 JacobianForm;
typedef Form_0::TestSpace FunctionSpace;

}

#endif
