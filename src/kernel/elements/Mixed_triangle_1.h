// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.3.5.

#ifndef __MIXED_TRIANGLE_1_H
#define __MIXED_TRIANGLE_1_H

#include <cmath>
#include <ufc.h>

/// This class defines the interface for a finite element.

class Mixed_triangle_1_finite_element_0_0: public ufc::finite_element
{
public:

  /// Constructor
  Mixed_triangle_1_finite_element_0_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Mixed_triangle_1_finite_element_0_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Discontinuous Lagrange finite element of degree 1 on a triangle";
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
    
    const static unsigned int dof = i;
    
    // Table(s) of coefficients
    const static double coefficients0[3][3] = \
    {{0.471404520791, -0.2886751345948, -0.1666666666667},
    {0.471404520791, 0.2886751345948, -0.1666666666667},
    {0.471404520791, 0, 0.3333333333333}};
    
    // Generate scalings
    const double scalings_y_0 = 1.0;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5 * y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1.0;
    const double psitilde_a_1 = 1*x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1.0;
    const double psitilde_bs_0_1 = 0.5 + 1.5*y;
    const double psitilde_bs_1_0 = 1.0;
    
    // Compute basisvalues
    const double basisvalues[3] = \
    {psitilde_a_0*scalings_y_0*psitilde_bs_0_0*0.7071067811865,
    psitilde_a_1*scalings_y_1*psitilde_bs_1_0*1.732050807569,
    psitilde_a_0*scalings_y_0*psitilde_bs_0_1*1};
    
    // Compute value(s)
    *values = 0.0;
    for (unsigned int j = 0; j < 3; j++)
      *values += coefficients0[dof][j]*basisvalues[j];
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    // Not implemented (only for Lagrange elements
    return 0;
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
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
    return new Mixed_triangle_1_finite_element_0_0();
  }

};

/// This class defines the interface for a finite element.

class Mixed_triangle_1_finite_element_0_1: public ufc::finite_element
{
public:

  /// Constructor
  Mixed_triangle_1_finite_element_0_1() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Mixed_triangle_1_finite_element_0_1()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Discontinuous Lagrange finite element of degree 1 on a triangle";
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
    
    const static unsigned int dof = i;
    
    // Table(s) of coefficients
    const static double coefficients0[3][3] = \
    {{0.471404520791, -0.2886751345948, -0.1666666666667},
    {0.471404520791, 0.2886751345948, -0.1666666666667},
    {0.471404520791, 0, 0.3333333333333}};
    
    // Generate scalings
    const double scalings_y_0 = 1.0;
    const double scalings_y_1 = scalings_y_0*(0.5 - 0.5 * y);
    
    // Compute psitilde_a
    const double psitilde_a_0 = 1.0;
    const double psitilde_a_1 = 1*x;
    
    // Compute psitilde_bs
    const double psitilde_bs_0_0 = 1.0;
    const double psitilde_bs_0_1 = 0.5 + 1.5*y;
    const double psitilde_bs_1_0 = 1.0;
    
    // Compute basisvalues
    const double basisvalues[3] = \
    {psitilde_a_0*scalings_y_0*psitilde_bs_0_0*0.7071067811865,
    psitilde_a_1*scalings_y_1*psitilde_bs_1_0*1.732050807569,
    psitilde_a_0*scalings_y_0*psitilde_bs_0_1*1};
    
    // Compute value(s)
    *values = 0.0;
    for (unsigned int j = 0; j < 3; j++)
      *values += coefficients0[dof][j]*basisvalues[j];
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    // Not implemented (only for Lagrange elements
    return 0;
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell& c) const
  {
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
    return new Mixed_triangle_1_finite_element_0_1();
  }

};

/// This class defines the interface for a finite element.

class Mixed_triangle_1_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  Mixed_triangle_1_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Mixed_triangle_1_finite_element_0()
  {
    // Do nothing
  }

  /// Return a string identifying the finite element
  virtual const char* signature() const
  {
    return "Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a triangle, Discontinuous Lagrange finite element of degree 1 on a triangle]";
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
    
    for (unsigned int element = 0; element < 2; element++)
    // Switch for each of the basis elements
    {
      switch ( element )
      {
        case 0:
        {
          if (0 <= i and i <= 2)
          {
            // Compute local degree of freedom
            const static unsigned int dof = i;
    
            // Table(s) of coefficients
            const static double coefficients0[3][3] =         \
            {{0.471404520791, -0.2886751345948, -0.1666666666667},
            {0.471404520791, 0.2886751345948, -0.1666666666667},
            {0.471404520791, 0, 0.3333333333333}};
    
            // Generate scalings
            const double scalings_y_0 = 1.0;
            const double scalings_y_1 = scalings_y_0*(0.5 - 0.5 * y);
    
            // Compute psitilde_a
            const double psitilde_a_0 = 1.0;
            const double psitilde_a_1 = 1*x;
    
            // Compute psitilde_bs
            const double psitilde_bs_0_0 = 1.0;
            const double psitilde_bs_0_1 = 0.5 + 1.5*y;
            const double psitilde_bs_1_0 = 1.0;
    
            // Compute basisvalues
            const double basisvalues[3] =         \
            {psitilde_a_0*scalings_y_0*psitilde_bs_0_0*0.7071067811865,
            psitilde_a_1*scalings_y_1*psitilde_bs_1_0*1.732050807569,
            psitilde_a_0*scalings_y_0*psitilde_bs_0_1*1};
    
            // Compute value(s)
            values[0] = 0.0;
            for (unsigned int j = 0; j < 3; j++)
            {
              values[0] += coefficients0[dof][j]*basisvalues[j];
            }
          }
          else
          {
            values[0] = 0.0;
          }
          break;
        }
        case 1:
        {
          if (3 <= i and i <= 5)
          {
            // Compute local degree of freedom
            const static unsigned int dof = i - 3;
    
            // Table(s) of coefficients
            const static double coefficients0[3][3] =         \
            {{0.471404520791, -0.2886751345948, -0.1666666666667},
            {0.471404520791, 0.2886751345948, -0.1666666666667},
            {0.471404520791, 0, 0.3333333333333}};
    
            // Generate scalings
            const double scalings_y_0 = 1.0;
            const double scalings_y_1 = scalings_y_0*(0.5 - 0.5 * y);
    
            // Compute psitilde_a
            const double psitilde_a_0 = 1.0;
            const double psitilde_a_1 = 1*x;
    
            // Compute psitilde_bs
            const double psitilde_bs_0_0 = 1.0;
            const double psitilde_bs_0_1 = 0.5 + 1.5*y;
            const double psitilde_bs_1_0 = 1.0;
    
            // Compute basisvalues
            const double basisvalues[3] =         \
            {psitilde_a_0*scalings_y_0*psitilde_bs_0_0*0.7071067811865,
            psitilde_a_1*scalings_y_1*psitilde_bs_1_0*1.732050807569,
            psitilde_a_0*scalings_y_0*psitilde_bs_0_1*1};
    
            // Compute value(s)
            values[1] = 0.0;
            for (unsigned int j = 0; j < 3; j++)
            {
              values[1] += coefficients0[dof][j]*basisvalues[j];
            }
          }
          else
          {
            values[1] = 0.0;
          }
          break;
        }
      }
    }
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    double values[2];
    double coordinates[2];
    
    // Nodal coordinates on reference cell
    static double X[6][2] = {{0, 0}, {1, 0}, {0, 1}, {0, 0}, {1, 0}, {0, 1}};
    
    // Components for each dof
    static unsigned int components[6] = {0, 0, 0, 1, 1, 1};
    
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
    vertex_values[0] = dof_values[0];
    vertex_values[1] = dof_values[1];
    vertex_values[2] = dof_values[2];
    vertex_values[3] = dof_values[3];
    vertex_values[4] = dof_values[4];
    vertex_values[5] = dof_values[5];
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 2;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    switch ( i )
    {
    case 0:
      return new Mixed_triangle_1_finite_element_0_0();
      break;
    case 1:
      return new Mixed_triangle_1_finite_element_0_1();
      break;
    }
    return 0;
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Mixed_triangle_1_dof_map_0_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Mixed_triangle_1_dof_map_0_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Mixed_triangle_1_dof_map_0_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Discontinuous Lagrange finite element of degree 1 on a triangle";
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
      return true;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = 3*m.num_entities[2];
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
    return 3;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 0;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
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
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::mesh& m,
                                    const ufc::cell& c) const
  {
    // Not implemented
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    // Not implemented
    return 0;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    // Not implemented
    return 0;
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Mixed_triangle_1_dof_map_0_1: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Mixed_triangle_1_dof_map_0_1() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Mixed_triangle_1_dof_map_0_1()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Discontinuous Lagrange finite element of degree 1 on a triangle";
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
      return true;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = 3*m.num_entities[2];
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
    return 3;
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 0;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
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
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::mesh& m,
                                    const ufc::cell& c) const
  {
    // Not implemented
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    // Not implemented
    return 0;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    // Not implemented
    return 0;
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Mixed_triangle_1_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Mixed_triangle_1_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Mixed_triangle_1_dof_map_0()
  {
    // Do nothing
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
    return "FFC dof map for Mixed finite element: [Discontinuous Lagrange finite element of degree 1 on a triangle, Discontinuous Lagrange finite element of degree 1 on a triangle]";
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
      return true;
      break;
    }
    return false;
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
    __global_dimension = 6*m.num_entities[2];
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

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
    return 0;
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
    dofs[0] = 3*c.entity_indices[2][0];
    dofs[1] = 3*c.entity_indices[2][0] + 1;
    dofs[2] = 3*c.entity_indices[2][0] + 2;
    unsigned int offset = 3*m.num_entities[2];
    dofs[3] = offset + 3*c.entity_indices[2][0];
    dofs[4] = offset + 3*c.entity_indices[2][0] + 1;
    dofs[5] = offset + 3*c.entity_indices[2][0] + 2;
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
      
      break;
    case 1:
      
      break;
    case 2:
      
      break;
    }
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::mesh& m,
                                    const ufc::cell& c) const
  {
    // Not implemented
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
    // Not implemented
    return 0;
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
    // Not implemented
    return 0;
  }

};

#endif
