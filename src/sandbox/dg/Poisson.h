// This code conforms with the UFC specification version 1.0
// and was automatically generated by FFC version 0.3.5.

#ifndef __POISSON_H
#define __POISSON_H

#include <ufc.h>

/// This class defines the interface for a finite element.

class Poisson_finite_element_0: public ufc::finite_element
{
public:

  /// Constructor
  Poisson_finite_element_0() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Poisson_finite_element_0()
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
    // Not implemented
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    // Not implemented
    return 0.0;
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell & c) const
  {
    // Not implemented
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new Poisson_finite_element_0();
  }

};

/// This class defines the interface for a finite element.

class Poisson_finite_element_1: public ufc::finite_element
{
public:

  /// Constructor
  Poisson_finite_element_1() : ufc::finite_element()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Poisson_finite_element_1()
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
    // Not implemented
  }

  /// Evaluate linear functional for dof i on the function f
  virtual double evaluate_dof(unsigned int i,
                              const ufc::function& f,
                              const ufc::cell& c) const
  {
    // Not implemented
    return 0.0;
  }

  /// Interpolate vertex values from dof values
  virtual void interpolate_vertex_values(double* vertex_values,
                                         const double* dof_values,
                                         const ufc::cell & c) const
  {
    // Not implemented
  }

  /// Return the number of sub elements (for a mixed element)
  virtual unsigned int num_sub_elements() const
  {
    return 1;
  }

  /// Create a new finite element for sub element i (for a mixed element)
  virtual ufc::finite_element* create_sub_element(unsigned int i) const
  {
    return new Poisson_finite_element_1();
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Poisson_dof_map_0: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Poisson_dof_map_0() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Poisson_dof_map_0()
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

  /// Return the number of dofs on a facets of a cell
  virtual unsigned int num_facet_dofs() const
  {
    // Not implemented
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

  /// Tabulate the local-to-global mapping of dofs on a facet of a cell
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   const ufc::mesh& m,
                                   const ufc::cell& c,
                                   unsigned int facet) const
  {
    // Not implemented
  }

};

/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class Poisson_dof_map_1: public ufc::dof_map
{
private:

  unsigned int __global_dimension;

public:

  /// Constructor
  Poisson_dof_map_1() : ufc::dof_map()
  {
    __global_dimension = 0;
  }

  /// Destructor
  virtual ~Poisson_dof_map_1()
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

  /// Return the number of dofs on a facets of a cell
  virtual unsigned int num_facet_dofs() const
  {
    // Not implemented
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

  /// Tabulate the local-to-global mapping of dofs on a facet of a cell
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   const ufc::mesh& m,
                                   const ufc::cell& c,
                                   unsigned int facet) const
  {
    // Not implemented
  }

};

/// This class defines the interface for the tabulation of the
/// interior facet tensor corresponding to the local contribution to
/// a form from the integral over an interior facet.

class Poisson_interior_facet_integral_0: public ufc::interior_facet_integral
{
public:

  /// Constructor
  Poisson_interior_facet_integral_0() : ufc::interior_facet_integral()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Poisson_interior_facet_integral_0()
  {
    // Do nothing
  }

  /// Tabulate the tensor for the contribution from a local interior facet
  virtual void tabulate_tensor(double* A,
                               const double * const * w,
                               const ufc::cell& c0,
                               const ufc::cell& c1,
                               unsigned int facet0,
                               unsigned int facet1) const
  {
    // Extract coordinates
    const double * const * x0 = c0.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J0_00 = x0[1][0] - x0[0][0];
    const double J0_01 = x0[2][0] - x0[0][0];
    const double J0_10 = x0[1][1] - x0[0][1];
    const double J0_11 = x0[2][1] - x0[0][1];
      
    // Compute determinant of Jacobian
    double detJ0 = J0_00*J0_11 - J0_01*J0_10;
      
    // Compute inverse of Jacobian
    // const double Jinv0_00 =  J0_11 / detJ0;
    // const double Jinv0_01 = -J0_01 / detJ0;
    // const double Jinv0_10 = -J0_10 / detJ0;
    // const double Jinv0_11 =  J0_00 / detJ0;
    
    // Take absolute value of determinant
    detJ0 = std::abs(detJ0);
    
    // Extract coordinates
    const double * const * x1 = c1.coordinates;
    
    // Compute Jacobian of affine map from reference cell
    const double J1_00 = x1[1][0] - x1[0][0];
    const double J1_01 = x1[2][0] - x1[0][0];
    const double J1_10 = x1[1][1] - x1[0][1];
    const double J1_11 = x1[2][1] - x1[0][1];
      
    // Compute determinant of Jacobian
    double detJ1 = J1_00*J1_11 - J1_01*J1_10;
      
    // Compute inverse of Jacobian
    // const double Jinv1_00 =  J1_11 / detJ1;
    // const double Jinv1_01 = -J1_01 / detJ1;
    // const double Jinv1_10 = -J1_10 / detJ1;
    // const double Jinv1_11 =  J1_00 / detJ1;
    
    // Take absolute value of determinant
    detJ1 = std::abs(detJ1);
    
    // Vertices on edges
    static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    
    // Get vertices
    const unsigned int v0 = edge_vertices[facet0][0];
    const unsigned int v1 = edge_vertices[facet0][1];
    
    // Compute scale factor (length of edge scaled by length of reference interval)
    const double dx0 = x0[v1][0] - x0[v0][0];
    const double dx1 = x0[v1][1] - x0[v0][1];
    const double det = std::sqrt(dx0*dx0 + dx1*dx1);
    
    // Compute geometry tensors
    const double G0_ = det;
    
    // Compute element tensor for all facet-facet combinations
    switch ( facet0 )
    {
    case 0:
      switch ( facet1 )
      {
      case 0:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 0.000000000000000e+00;
        A[4] = 0.000000000000000e+00;
        A[5] = 0.000000000000000e+00;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 0.000000000000000e+00;
        A[10] = 3.333333333333331e-01*G0_;
        A[11] = 1.666666666666666e-01*G0_;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 0.000000000000000e+00;
        A[16] = 1.666666666666665e-01*G0_;
        A[17] = 3.333333333333330e-01*G0_;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      case 1:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 0.000000000000000e+00;
        A[4] = 0.000000000000000e+00;
        A[5] = 0.000000000000000e+00;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 1.666666666666666e-01*G0_;
        A[10] = 0.000000000000000e+00;
        A[11] = 3.333333333333330e-01*G0_;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 3.333333333333330e-01*G0_;
        A[16] = 0.000000000000000e+00;
        A[17] = 1.666666666666665e-01*G0_;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      case 2:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 0.000000000000000e+00;
        A[4] = 0.000000000000000e+00;
        A[5] = 0.000000000000000e+00;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 3.333333333333330e-01*G0_;
        A[10] = 1.666666666666666e-01*G0_;
        A[11] = 0.000000000000000e+00;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 1.666666666666665e-01*G0_;
        A[16] = 3.333333333333330e-01*G0_;
        A[17] = 0.000000000000000e+00;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      }
      
      break;
    case 1:
      switch ( facet1 )
      {
      case 0:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 0.000000000000000e+00;
        A[4] = 1.666666666666665e-01*G0_;
        A[5] = 3.333333333333330e-01*G0_;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 0.000000000000000e+00;
        A[10] = 0.000000000000000e+00;
        A[11] = 0.000000000000000e+00;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 0.000000000000000e+00;
        A[16] = 3.333333333333330e-01*G0_;
        A[17] = 1.666666666666665e-01*G0_;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      case 1:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 3.333333333333330e-01*G0_;
        A[4] = 0.000000000000000e+00;
        A[5] = 1.666666666666665e-01*G0_;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 0.000000000000000e+00;
        A[10] = 0.000000000000000e+00;
        A[11] = 0.000000000000000e+00;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 1.666666666666665e-01*G0_;
        A[16] = 0.000000000000000e+00;
        A[17] = 3.333333333333330e-01*G0_;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      case 2:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 1.666666666666665e-01*G0_;
        A[4] = 3.333333333333330e-01*G0_;
        A[5] = 0.000000000000000e+00;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 0.000000000000000e+00;
        A[10] = 0.000000000000000e+00;
        A[11] = 0.000000000000000e+00;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 3.333333333333330e-01*G0_;
        A[16] = 1.666666666666665e-01*G0_;
        A[17] = 0.000000000000000e+00;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      }
      
      break;
    case 2:
      switch ( facet1 )
      {
      case 0:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 0.000000000000000e+00;
        A[4] = 3.333333333333330e-01*G0_;
        A[5] = 1.666666666666665e-01*G0_;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 0.000000000000000e+00;
        A[10] = 1.666666666666666e-01*G0_;
        A[11] = 3.333333333333330e-01*G0_;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 0.000000000000000e+00;
        A[16] = 0.000000000000000e+00;
        A[17] = 0.000000000000000e+00;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      case 1:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 1.666666666666665e-01*G0_;
        A[4] = 0.000000000000000e+00;
        A[5] = 3.333333333333330e-01*G0_;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 3.333333333333330e-01*G0_;
        A[10] = 0.000000000000000e+00;
        A[11] = 1.666666666666665e-01*G0_;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 0.000000000000000e+00;
        A[16] = 0.000000000000000e+00;
        A[17] = 0.000000000000000e+00;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      case 2:
        A[0] = 0.000000000000000e+00;
        A[1] = 0.000000000000000e+00;
        A[2] = 0.000000000000000e+00;
        A[3] = 3.333333333333330e-01*G0_;
        A[4] = 1.666666666666665e-01*G0_;
        A[5] = 0.000000000000000e+00;
        A[6] = 0.000000000000000e+00;
        A[7] = 0.000000000000000e+00;
        A[8] = 0.000000000000000e+00;
        A[9] = 1.666666666666665e-01*G0_;
        A[10] = 3.333333333333331e-01*G0_;
        A[11] = 0.000000000000000e+00;
        A[12] = 0.000000000000000e+00;
        A[13] = 0.000000000000000e+00;
        A[14] = 0.000000000000000e+00;
        A[15] = 0.000000000000000e+00;
        A[16] = 0.000000000000000e+00;
        A[17] = 0.000000000000000e+00;
        A[18] = 0.000000000000000e+00;
        A[19] = 0.000000000000000e+00;
        A[20] = 0.000000000000000e+00;
        A[21] = 0.000000000000000e+00;
        A[22] = 0.000000000000000e+00;
        A[23] = 0.000000000000000e+00;
        A[24] = 0.000000000000000e+00;
        A[25] = 0.000000000000000e+00;
        A[26] = 0.000000000000000e+00;
        A[27] = 0.000000000000000e+00;
        A[28] = 0.000000000000000e+00;
        A[29] = 0.000000000000000e+00;
        A[30] = 0.000000000000000e+00;
        A[31] = 0.000000000000000e+00;
        A[32] = 0.000000000000000e+00;
        A[33] = 0.000000000000000e+00;
        A[34] = 0.000000000000000e+00;
        A[35] = 0.000000000000000e+00;
        break;
      }
      
      break;
    }
    
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

class Poisson: public ufc::form
{
public:

  /// Constructor
  Poisson() : ufc::form()
  {
    // Do nothing
  }

  /// Destructor
  virtual ~Poisson()
  {
    // Do nothing
  }

  /// Return a string identifying the form
  virtual const char* signature() const
  {
    return "|det(F)|^(1) | vi0(+)*vi1(-)*dS(0)";
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
    return 1;
  }
    
  /// Create a new finite element for argument function i
  virtual ufc::finite_element* create_finite_element(unsigned int i) const
  {
    switch ( i )
    {
    case 0:
      return new Poisson_finite_element_0();
      break;
    case 1:
      return new Poisson_finite_element_1();
      break;
    }
    return 0;
  }
  
  /// Create a new dof map for argument function i
  virtual ufc::dof_map* create_dof_map(unsigned int i) const
  {
    switch ( i )
    {
    case 0:
      return new Poisson_dof_map_0();
      break;
    case 1:
      return new Poisson_dof_map_1();
      break;
    }
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
    return new Poisson_interior_facet_integral_0();
  }

};

#endif
