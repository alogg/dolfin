// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.4.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __ENERGYNORM_H
#define __ENERGYNORM_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>
#include <dolfin/BilinearForm.h>
#include <dolfin/LinearForm.h>
#include <dolfin/Functional.h>
#include <dolfin/FEM.h>

namespace dolfin { namespace EnergyNorm {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class Functional : public dolfin::Functional
{
public:

  class FunctionElement_0;

  Functional();
  
  real operator() (Function& w0, Mesh& mesh);

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

};

class Functional::FunctionElement_0 : public dolfin::FiniteElement
{
public:

  FunctionElement_0() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_0()
  {
    if ( tensordims ) delete [] tensordims;
    if ( subelements )
    {
      for (unsigned int i = 0; i < elementdim(); i++)
        delete subelements[i];
      delete [] subelements;
    }
  }

  inline unsigned int spacedim() const
  {
    return 6;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    dolfin_error("Element is scalar.");
    return 0;
  }

  inline unsigned int elementdim() const
  {
    return 1;
  }

  inline unsigned int rank() const
  {
    return 0;
  }

  void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
  {
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int offset = mesh.topology().size(0);
    nodes[3] = offset + cell.entities(1)[0];
    nodes[4] = offset + cell.entities(1)[1];
    nodes[5] = offset + cell.entities(1)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
    points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
    points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
  }

  const FiniteElement& operator[] (unsigned int i) const
  {
    return *this;
  }

  FiniteElement& operator[] (unsigned int i)
  {
    return *this;
  }

  FiniteElementSpec spec() const
  {
    FiniteElementSpec s("Lagrange", "triangle", 2);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

Functional::Functional() : dolfin::Functional(1)
{
}

real Functional::operator() (Function& w0, Mesh& mesh)
{
  // Initialize functions
  initFunction(0, w0, new FunctionElement_0());

  // Assemble value of functional
  return FEM::assemble(*this, mesh);
}

// Contribution from the interior
bool Functional::interior_contribution() const { return true; }

void Functional::eval(real block[], const AffineMap& map) const
{
  // Compute coefficients
  const real c0_0 = c[0][0];
  const real c0_1 = c[0][1];
  const real c0_2 = c[0][2];
  const real c0_3 = c[0][3];
  const real c0_4 = c[0][4];
  const real c0_5 = c[0][5];

  // Compute geometry tensors
  const real G0_0_0 = map.det*c0_0*c0_0;
  const real G0_0_1 = map.det*c0_0*c0_1;
  const real G0_0_2 = map.det*c0_0*c0_2;
  const real G0_0_3 = map.det*c0_0*c0_3;
  const real G0_1_0 = map.det*c0_1*c0_0;
  const real G0_1_1 = map.det*c0_1*c0_1;
  const real G0_1_2 = map.det*c0_1*c0_2;
  const real G0_1_4 = map.det*c0_1*c0_4;
  const real G0_2_0 = map.det*c0_2*c0_0;
  const real G0_2_1 = map.det*c0_2*c0_1;
  const real G0_2_2 = map.det*c0_2*c0_2;
  const real G0_2_5 = map.det*c0_2*c0_5;
  const real G0_3_0 = map.det*c0_3*c0_0;
  const real G0_3_3 = map.det*c0_3*c0_3;
  const real G0_3_4 = map.det*c0_3*c0_4;
  const real G0_3_5 = map.det*c0_3*c0_5;
  const real G0_4_1 = map.det*c0_4*c0_1;
  const real G0_4_3 = map.det*c0_4*c0_3;
  const real G0_4_4 = map.det*c0_4*c0_4;
  const real G0_4_5 = map.det*c0_4*c0_5;
  const real G0_5_2 = map.det*c0_5*c0_2;
  const real G0_5_3 = map.det*c0_5*c0_3;
  const real G0_5_4 = map.det*c0_5*c0_4;
  const real G0_5_5 = map.det*c0_5*c0_5;
  const real G1_0_0_0_0 = map.det*c0_0*c0_0*map.g00*map.g00 + map.det*c0_0*c0_0*map.g01*map.g01;
  const real G1_0_0_0_1 = map.det*c0_0*c0_0*map.g00*map.g10 + map.det*c0_0*c0_0*map.g01*map.g11;
  const real G1_0_0_1_0 = map.det*c0_0*c0_1*map.g00*map.g00 + map.det*c0_0*c0_1*map.g01*map.g01;
  const real G1_0_0_2_1 = map.det*c0_0*c0_2*map.g00*map.g10 + map.det*c0_0*c0_2*map.g01*map.g11;
  const real G1_0_0_4_1 = map.det*c0_0*c0_4*map.g00*map.g10 + map.det*c0_0*c0_4*map.g01*map.g11;
  const real G1_0_0_5_0 = map.det*c0_0*c0_5*map.g00*map.g00 + map.det*c0_0*c0_5*map.g01*map.g01;
  const real G1_0_1_0_0 = map.det*c0_0*c0_0*map.g10*map.g00 + map.det*c0_0*c0_0*map.g11*map.g01;
  const real G1_0_1_0_1 = map.det*c0_0*c0_0*map.g10*map.g10 + map.det*c0_0*c0_0*map.g11*map.g11;
  const real G1_0_1_1_0 = map.det*c0_0*c0_1*map.g10*map.g00 + map.det*c0_0*c0_1*map.g11*map.g01;
  const real G1_0_1_2_1 = map.det*c0_0*c0_2*map.g10*map.g10 + map.det*c0_0*c0_2*map.g11*map.g11;
  const real G1_0_1_4_1 = map.det*c0_0*c0_4*map.g10*map.g10 + map.det*c0_0*c0_4*map.g11*map.g11;
  const real G1_0_1_5_0 = map.det*c0_0*c0_5*map.g10*map.g00 + map.det*c0_0*c0_5*map.g11*map.g01;
  const real G1_1_0_0_0 = map.det*c0_1*c0_0*map.g00*map.g00 + map.det*c0_1*c0_0*map.g01*map.g01;
  const real G1_1_0_0_1 = map.det*c0_1*c0_0*map.g00*map.g10 + map.det*c0_1*c0_0*map.g01*map.g11;
  const real G1_1_0_1_0 = map.det*c0_1*c0_1*map.g00*map.g00 + map.det*c0_1*c0_1*map.g01*map.g01;
  const real G1_1_0_2_1 = map.det*c0_1*c0_2*map.g00*map.g10 + map.det*c0_1*c0_2*map.g01*map.g11;
  const real G1_1_0_3_1 = map.det*c0_1*c0_3*map.g00*map.g10 + map.det*c0_1*c0_3*map.g01*map.g11;
  const real G1_1_0_5_0 = map.det*c0_1*c0_5*map.g00*map.g00 + map.det*c0_1*c0_5*map.g01*map.g01;
  const real G1_1_0_5_1 = map.det*c0_1*c0_5*map.g00*map.g10 + map.det*c0_1*c0_5*map.g01*map.g11;
  const real G1_2_1_0_0 = map.det*c0_2*c0_0*map.g10*map.g00 + map.det*c0_2*c0_0*map.g11*map.g01;
  const real G1_2_1_0_1 = map.det*c0_2*c0_0*map.g10*map.g10 + map.det*c0_2*c0_0*map.g11*map.g11;
  const real G1_2_1_1_0 = map.det*c0_2*c0_1*map.g10*map.g00 + map.det*c0_2*c0_1*map.g11*map.g01;
  const real G1_2_1_2_1 = map.det*c0_2*c0_2*map.g10*map.g10 + map.det*c0_2*c0_2*map.g11*map.g11;
  const real G1_2_1_3_0 = map.det*c0_2*c0_3*map.g10*map.g00 + map.det*c0_2*c0_3*map.g11*map.g01;
  const real G1_2_1_4_0 = map.det*c0_2*c0_4*map.g10*map.g00 + map.det*c0_2*c0_4*map.g11*map.g01;
  const real G1_2_1_4_1 = map.det*c0_2*c0_4*map.g10*map.g10 + map.det*c0_2*c0_4*map.g11*map.g11;
  const real G1_3_0_2_1 = map.det*c0_3*c0_2*map.g00*map.g10 + map.det*c0_3*c0_2*map.g01*map.g11;
  const real G1_3_0_3_0 = map.det*c0_3*c0_3*map.g00*map.g00 + map.det*c0_3*c0_3*map.g01*map.g01;
  const real G1_3_0_3_1 = map.det*c0_3*c0_3*map.g00*map.g10 + map.det*c0_3*c0_3*map.g01*map.g11;
  const real G1_3_0_4_0 = map.det*c0_3*c0_4*map.g00*map.g00 + map.det*c0_3*c0_4*map.g01*map.g01;
  const real G1_3_0_4_1 = map.det*c0_3*c0_4*map.g00*map.g10 + map.det*c0_3*c0_4*map.g01*map.g11;
  const real G1_3_0_5_1 = map.det*c0_3*c0_5*map.g00*map.g10 + map.det*c0_3*c0_5*map.g01*map.g11;
  const real G1_3_1_1_0 = map.det*c0_3*c0_1*map.g10*map.g00 + map.det*c0_3*c0_1*map.g11*map.g01;
  const real G1_3_1_3_0 = map.det*c0_3*c0_3*map.g10*map.g00 + map.det*c0_3*c0_3*map.g11*map.g01;
  const real G1_3_1_3_1 = map.det*c0_3*c0_3*map.g10*map.g10 + map.det*c0_3*c0_3*map.g11*map.g11;
  const real G1_3_1_4_0 = map.det*c0_3*c0_4*map.g10*map.g00 + map.det*c0_3*c0_4*map.g11*map.g01;
  const real G1_3_1_5_0 = map.det*c0_3*c0_5*map.g10*map.g00 + map.det*c0_3*c0_5*map.g11*map.g01;
  const real G1_3_1_5_1 = map.det*c0_3*c0_5*map.g10*map.g10 + map.det*c0_3*c0_5*map.g11*map.g11;
  const real G1_4_0_2_1 = map.det*c0_4*c0_2*map.g00*map.g10 + map.det*c0_4*c0_2*map.g01*map.g11;
  const real G1_4_0_3_0 = map.det*c0_4*c0_3*map.g00*map.g00 + map.det*c0_4*c0_3*map.g01*map.g01;
  const real G1_4_0_3_1 = map.det*c0_4*c0_3*map.g00*map.g10 + map.det*c0_4*c0_3*map.g01*map.g11;
  const real G1_4_0_4_0 = map.det*c0_4*c0_4*map.g00*map.g00 + map.det*c0_4*c0_4*map.g01*map.g01;
  const real G1_4_0_4_1 = map.det*c0_4*c0_4*map.g00*map.g10 + map.det*c0_4*c0_4*map.g01*map.g11;
  const real G1_4_0_5_1 = map.det*c0_4*c0_5*map.g00*map.g10 + map.det*c0_4*c0_5*map.g01*map.g11;
  const real G1_4_1_0_0 = map.det*c0_4*c0_0*map.g10*map.g00 + map.det*c0_4*c0_0*map.g11*map.g01;
  const real G1_4_1_0_1 = map.det*c0_4*c0_0*map.g10*map.g10 + map.det*c0_4*c0_0*map.g11*map.g11;
  const real G1_4_1_2_1 = map.det*c0_4*c0_2*map.g10*map.g10 + map.det*c0_4*c0_2*map.g11*map.g11;
  const real G1_4_1_3_0 = map.det*c0_4*c0_3*map.g10*map.g00 + map.det*c0_4*c0_3*map.g11*map.g01;
  const real G1_4_1_4_0 = map.det*c0_4*c0_4*map.g10*map.g00 + map.det*c0_4*c0_4*map.g11*map.g01;
  const real G1_4_1_4_1 = map.det*c0_4*c0_4*map.g10*map.g10 + map.det*c0_4*c0_4*map.g11*map.g11;
  const real G1_4_1_5_0 = map.det*c0_4*c0_5*map.g10*map.g00 + map.det*c0_4*c0_5*map.g11*map.g01;
  const real G1_5_0_0_0 = map.det*c0_5*c0_0*map.g00*map.g00 + map.det*c0_5*c0_0*map.g01*map.g01;
  const real G1_5_0_0_1 = map.det*c0_5*c0_0*map.g00*map.g10 + map.det*c0_5*c0_0*map.g01*map.g11;
  const real G1_5_0_1_0 = map.det*c0_5*c0_1*map.g00*map.g00 + map.det*c0_5*c0_1*map.g01*map.g01;
  const real G1_5_0_3_1 = map.det*c0_5*c0_3*map.g00*map.g10 + map.det*c0_5*c0_3*map.g01*map.g11;
  const real G1_5_0_4_1 = map.det*c0_5*c0_4*map.g00*map.g10 + map.det*c0_5*c0_4*map.g01*map.g11;
  const real G1_5_0_5_0 = map.det*c0_5*c0_5*map.g00*map.g00 + map.det*c0_5*c0_5*map.g01*map.g01;
  const real G1_5_0_5_1 = map.det*c0_5*c0_5*map.g00*map.g10 + map.det*c0_5*c0_5*map.g01*map.g11;
  const real G1_5_1_1_0 = map.det*c0_5*c0_1*map.g10*map.g00 + map.det*c0_5*c0_1*map.g11*map.g01;
  const real G1_5_1_3_0 = map.det*c0_5*c0_3*map.g10*map.g00 + map.det*c0_5*c0_3*map.g11*map.g01;
  const real G1_5_1_3_1 = map.det*c0_5*c0_3*map.g10*map.g10 + map.det*c0_5*c0_3*map.g11*map.g11;
  const real G1_5_1_4_0 = map.det*c0_5*c0_4*map.g10*map.g00 + map.det*c0_5*c0_4*map.g11*map.g01;
  const real G1_5_1_5_0 = map.det*c0_5*c0_5*map.g10*map.g00 + map.det*c0_5*c0_5*map.g11*map.g01;
  const real G1_5_1_5_1 = map.det*c0_5*c0_5*map.g10*map.g10 + map.det*c0_5*c0_5*map.g11*map.g11;

  // Compute element tensor
  block[0] = 1.666666666666665e-02*G0_0_0 - 2.777777777777774e-03*G0_0_1 - 2.777777777777775e-03*G0_0_2 - 1.111111111111110e-02*G0_0_3 - 2.777777777777774e-03*G0_1_0 + 1.666666666666665e-02*G0_1_1 - 2.777777777777776e-03*G0_1_2 - 1.111111111111111e-02*G0_1_4 - 2.777777777777775e-03*G0_2_0 - 2.777777777777776e-03*G0_2_1 + 1.666666666666666e-02*G0_2_2 - 1.111111111111111e-02*G0_2_5 - 1.111111111111110e-02*G0_3_0 + 8.888888888888882e-02*G0_3_3 + 4.444444444444443e-02*G0_3_4 + 4.444444444444443e-02*G0_3_5 - 1.111111111111111e-02*G0_4_1 + 4.444444444444443e-02*G0_4_3 + 8.888888888888884e-02*G0_4_4 + 4.444444444444442e-02*G0_4_5 - 1.111111111111111e-02*G0_5_2 + 4.444444444444443e-02*G0_5_3 + 4.444444444444443e-02*G0_5_4 + 8.888888888888882e-02*G0_5_5 + 4.999999999999992e-01*G1_0_0_0_0 + 4.999999999999991e-01*G1_0_0_0_1 + 1.666666666666663e-01*G1_0_0_1_0 + 1.666666666666665e-01*G1_0_0_2_1 - 6.666666666666657e-01*G1_0_0_4_1 - 6.666666666666655e-01*G1_0_0_5_0 + 4.999999999999991e-01*G1_0_1_0_0 + 4.999999999999991e-01*G1_0_1_0_1 + 1.666666666666664e-01*G1_0_1_1_0 + 1.666666666666665e-01*G1_0_1_2_1 - 6.666666666666656e-01*G1_0_1_4_1 - 6.666666666666655e-01*G1_0_1_5_0 + 1.666666666666663e-01*G1_1_0_0_0 + 1.666666666666664e-01*G1_1_0_0_1 + 4.999999999999992e-01*G1_1_0_1_0 - 1.666666666666664e-01*G1_1_0_2_1 + 6.666666666666654e-01*G1_1_0_3_1 - 6.666666666666654e-01*G1_1_0_5_0 - 6.666666666666656e-01*G1_1_0_5_1 + 1.666666666666665e-01*G1_2_1_0_0 + 1.666666666666665e-01*G1_2_1_0_1 - 1.666666666666664e-01*G1_2_1_1_0 + 4.999999999999991e-01*G1_2_1_2_1 + 6.666666666666653e-01*G1_2_1_3_0 - 6.666666666666654e-01*G1_2_1_4_0 - 6.666666666666656e-01*G1_2_1_4_1 + 6.666666666666653e-01*G1_3_0_2_1 + 1.333333333333330e+00*G1_3_0_3_0 + 6.666666666666652e-01*G1_3_0_3_1 - 1.333333333333331e+00*G1_3_0_4_0 - 6.666666666666656e-01*G1_3_0_4_1 - 6.666666666666652e-01*G1_3_0_5_1 + 6.666666666666654e-01*G1_3_1_1_0 + 6.666666666666652e-01*G1_3_1_3_0 + 1.333333333333330e+00*G1_3_1_3_1 - 6.666666666666653e-01*G1_3_1_4_0 - 6.666666666666653e-01*G1_3_1_5_0 - 1.333333333333331e+00*G1_3_1_5_1 - 6.666666666666654e-01*G1_4_0_2_1 - 1.333333333333331e+00*G1_4_0_3_0 - 6.666666666666653e-01*G1_4_0_3_1 + 1.333333333333331e+00*G1_4_0_4_0 + 6.666666666666656e-01*G1_4_0_4_1 + 6.666666666666654e-01*G1_4_0_5_1 - 6.666666666666656e-01*G1_4_1_0_0 - 6.666666666666656e-01*G1_4_1_0_1 - 6.666666666666659e-01*G1_4_1_2_1 - 6.666666666666656e-01*G1_4_1_3_0 + 6.666666666666657e-01*G1_4_1_4_0 + 1.333333333333331e+00*G1_4_1_4_1 + 6.666666666666655e-01*G1_4_1_5_0 - 6.666666666666655e-01*G1_5_0_0_0 - 6.666666666666656e-01*G1_5_0_0_1 - 6.666666666666654e-01*G1_5_0_1_0 - 6.666666666666652e-01*G1_5_0_3_1 + 6.666666666666657e-01*G1_5_0_4_1 + 1.333333333333331e+00*G1_5_0_5_0 + 6.666666666666656e-01*G1_5_0_5_1 - 6.666666666666656e-01*G1_5_1_1_0 - 6.666666666666653e-01*G1_5_1_3_0 - 1.333333333333331e+00*G1_5_1_3_1 + 6.666666666666654e-01*G1_5_1_4_0 + 6.666666666666656e-01*G1_5_1_5_0 + 1.333333333333331e+00*G1_5_1_5_1;
}

// No contribution from the boundary
bool Functional::boundary_contribution() const { return false; }

void Functional::eval(real block[], const AffineMap& map, unsigned int facet) const {}

} }

#endif
