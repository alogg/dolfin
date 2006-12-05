// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __L2ERROR_H
#define __L2ERROR_H

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

namespace dolfin { namespace L2Error {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:

  class TestElement;

  class FunctionElement_0;

  class FunctionElement_1;

  LinearForm(Function& w0, Function& w1);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

};

class LinearForm::TestElement : public dolfin::FiniteElement
{
public:

  TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~TestElement()
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
    return 1;
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
    nodes[0] = cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(3.333333333333334e-01, 3.333333333333334e-01);
    components[0] = 0;
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
    FiniteElementSpec s("Discontinuous Lagrange", "triangle", 0);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_0 : public dolfin::FiniteElement
{
public:

  FunctionElement_0() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 2;

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
    return 12;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    dolfin_assert(i < 1);
    return tensordims[i];
  }

  inline unsigned int elementdim() const
  {
    return 1;
  }

  inline unsigned int rank() const
  {
    return 1;
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
    offset = offset + mesh.topology().size(1);
    nodes[6] = offset + cell.entities(0)[0];
    nodes[7] = offset + cell.entities(0)[1];
    nodes[8] = offset + cell.entities(0)[2];
    offset = offset + mesh.topology().size(0);
    nodes[9] = offset + cell.entities(1)[0];
    nodes[10] = offset + cell.entities(1)[1];
    nodes[11] = offset + cell.entities(1)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
    points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
    points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[8] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[9] = map(5.000000000000000e-01, 5.000000000000000e-01);
    points[10] = map(0.000000000000000e+00, 5.000000000000000e-01);
    points[11] = map(5.000000000000000e-01, 0.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
    components[6] = 1;
    components[7] = 1;
    components[8] = 1;
    components[9] = 1;
    components[10] = 1;
    components[11] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0) + mesh.topology().size(1);
    vertex_nodes[1] = offset + vertex;
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 2, 2);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_1 : public dolfin::FiniteElement
{
public:

  FunctionElement_1() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 2;

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_1()
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
    return 12;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    dolfin_assert(i < 1);
    return tensordims[i];
  }

  inline unsigned int elementdim() const
  {
    return 1;
  }

  inline unsigned int rank() const
  {
    return 1;
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
    offset = offset + mesh.topology().size(1);
    nodes[6] = offset + cell.entities(0)[0];
    nodes[7] = offset + cell.entities(0)[1];
    nodes[8] = offset + cell.entities(0)[2];
    offset = offset + mesh.topology().size(0);
    nodes[9] = offset + cell.entities(1)[0];
    nodes[10] = offset + cell.entities(1)[1];
    nodes[11] = offset + cell.entities(1)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
    points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
    points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[8] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[9] = map(5.000000000000000e-01, 5.000000000000000e-01);
    points[10] = map(0.000000000000000e+00, 5.000000000000000e-01);
    points[11] = map(5.000000000000000e-01, 0.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
    components[6] = 1;
    components[7] = 1;
    components[8] = 1;
    components[9] = 1;
    components[10] = 1;
    components[11] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0) + mesh.topology().size(1);
    vertex_nodes[1] = offset + vertex;
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 2, 2);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

LinearForm::LinearForm(Function& w0, Function& w1) : dolfin::LinearForm(2)
{
  // Create finite element for test space
  _test = new TestElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
  initFunction(1, w1, new FunctionElement_1());
}

// Contribution from the interior
bool LinearForm::interior_contribution() const { return true; }

void LinearForm::eval(real block[], const AffineMap& map) const
{
  // Compute coefficients
  const real c0_0 = c[0][0];
  const real c0_1 = c[0][1];
  const real c0_2 = c[0][2];
  const real c0_3 = c[0][3];
  const real c0_4 = c[0][4];
  const real c0_5 = c[0][5];
  const real c0_6 = c[0][6];
  const real c0_7 = c[0][7];
  const real c0_8 = c[0][8];
  const real c0_9 = c[0][9];
  const real c0_10 = c[0][10];
  const real c0_11 = c[0][11];
  const real c1_0 = c[1][0];
  const real c1_1 = c[1][1];
  const real c1_2 = c[1][2];
  const real c1_3 = c[1][3];
  const real c1_4 = c[1][4];
  const real c1_5 = c[1][5];
  const real c1_6 = c[1][6];
  const real c1_7 = c[1][7];
  const real c1_8 = c[1][8];
  const real c1_9 = c[1][9];
  const real c1_10 = c[1][10];
  const real c1_11 = c[1][11];

  // Compute geometry tensors
  const real G0_0_0 = map.det*c0_0*c0_0 + map.det*c1_0*c1_0;
  const real G0_0_1 = map.det*c0_0*c0_1 + map.det*c1_0*c1_1;
  const real G0_0_2 = map.det*c0_0*c0_2 + map.det*c1_0*c1_2;
  const real G0_0_3 = map.det*c0_0*c0_3 + map.det*c1_0*c1_3;
  const real G0_1_0 = map.det*c0_1*c0_0 + map.det*c1_1*c1_0;
  const real G0_1_1 = map.det*c0_1*c0_1 + map.det*c1_1*c1_1;
  const real G0_1_2 = map.det*c0_1*c0_2 + map.det*c1_1*c1_2;
  const real G0_1_4 = map.det*c0_1*c0_4 + map.det*c1_1*c1_4;
  const real G0_2_0 = map.det*c0_2*c0_0 + map.det*c1_2*c1_0;
  const real G0_2_1 = map.det*c0_2*c0_1 + map.det*c1_2*c1_1;
  const real G0_2_2 = map.det*c0_2*c0_2 + map.det*c1_2*c1_2;
  const real G0_2_5 = map.det*c0_2*c0_5 + map.det*c1_2*c1_5;
  const real G0_3_0 = map.det*c0_3*c0_0 + map.det*c1_3*c1_0;
  const real G0_3_3 = map.det*c0_3*c0_3 + map.det*c1_3*c1_3;
  const real G0_3_4 = map.det*c0_3*c0_4 + map.det*c1_3*c1_4;
  const real G0_3_5 = map.det*c0_3*c0_5 + map.det*c1_3*c1_5;
  const real G0_4_1 = map.det*c0_4*c0_1 + map.det*c1_4*c1_1;
  const real G0_4_3 = map.det*c0_4*c0_3 + map.det*c1_4*c1_3;
  const real G0_4_4 = map.det*c0_4*c0_4 + map.det*c1_4*c1_4;
  const real G0_4_5 = map.det*c0_4*c0_5 + map.det*c1_4*c1_5;
  const real G0_5_2 = map.det*c0_5*c0_2 + map.det*c1_5*c1_2;
  const real G0_5_3 = map.det*c0_5*c0_3 + map.det*c1_5*c1_3;
  const real G0_5_4 = map.det*c0_5*c0_4 + map.det*c1_5*c1_4;
  const real G0_5_5 = map.det*c0_5*c0_5 + map.det*c1_5*c1_5;
  const real G0_6_6 = map.det*c0_6*c0_6 + map.det*c1_6*c1_6;
  const real G0_6_7 = map.det*c0_6*c0_7 + map.det*c1_6*c1_7;
  const real G0_6_8 = map.det*c0_6*c0_8 + map.det*c1_6*c1_8;
  const real G0_6_9 = map.det*c0_6*c0_9 + map.det*c1_6*c1_9;
  const real G0_7_6 = map.det*c0_7*c0_6 + map.det*c1_7*c1_6;
  const real G0_7_7 = map.det*c0_7*c0_7 + map.det*c1_7*c1_7;
  const real G0_7_8 = map.det*c0_7*c0_8 + map.det*c1_7*c1_8;
  const real G0_7_10 = map.det*c0_7*c0_10 + map.det*c1_7*c1_10;
  const real G0_8_6 = map.det*c0_8*c0_6 + map.det*c1_8*c1_6;
  const real G0_8_7 = map.det*c0_8*c0_7 + map.det*c1_8*c1_7;
  const real G0_8_8 = map.det*c0_8*c0_8 + map.det*c1_8*c1_8;
  const real G0_8_11 = map.det*c0_8*c0_11 + map.det*c1_8*c1_11;
  const real G0_9_6 = map.det*c0_9*c0_6 + map.det*c1_9*c1_6;
  const real G0_9_9 = map.det*c0_9*c0_9 + map.det*c1_9*c1_9;
  const real G0_9_10 = map.det*c0_9*c0_10 + map.det*c1_9*c1_10;
  const real G0_9_11 = map.det*c0_9*c0_11 + map.det*c1_9*c1_11;
  const real G0_10_7 = map.det*c0_10*c0_7 + map.det*c1_10*c1_7;
  const real G0_10_9 = map.det*c0_10*c0_9 + map.det*c1_10*c1_9;
  const real G0_10_10 = map.det*c0_10*c0_10 + map.det*c1_10*c1_10;
  const real G0_10_11 = map.det*c0_10*c0_11 + map.det*c1_10*c1_11;
  const real G0_11_8 = map.det*c0_11*c0_8 + map.det*c1_11*c1_8;
  const real G0_11_9 = map.det*c0_11*c0_9 + map.det*c1_11*c1_9;
  const real G0_11_10 = map.det*c0_11*c0_10 + map.det*c1_11*c1_10;
  const real G0_11_11 = map.det*c0_11*c0_11 + map.det*c1_11*c1_11;
  const real G1_0_0 = map.det*c0_0*c1_0 + map.det*c1_0*c0_0;
  const real G1_0_1 = map.det*c0_0*c1_1 + map.det*c1_0*c0_1;
  const real G1_0_2 = map.det*c0_0*c1_2 + map.det*c1_0*c0_2;
  const real G1_0_3 = map.det*c0_0*c1_3 + map.det*c1_0*c0_3;
  const real G1_1_0 = map.det*c0_1*c1_0 + map.det*c1_1*c0_0;
  const real G1_1_1 = map.det*c0_1*c1_1 + map.det*c1_1*c0_1;
  const real G1_1_2 = map.det*c0_1*c1_2 + map.det*c1_1*c0_2;
  const real G1_1_4 = map.det*c0_1*c1_4 + map.det*c1_1*c0_4;
  const real G1_2_0 = map.det*c0_2*c1_0 + map.det*c1_2*c0_0;
  const real G1_2_1 = map.det*c0_2*c1_1 + map.det*c1_2*c0_1;
  const real G1_2_2 = map.det*c0_2*c1_2 + map.det*c1_2*c0_2;
  const real G1_2_5 = map.det*c0_2*c1_5 + map.det*c1_2*c0_5;
  const real G1_3_0 = map.det*c0_3*c1_0 + map.det*c1_3*c0_0;
  const real G1_3_3 = map.det*c0_3*c1_3 + map.det*c1_3*c0_3;
  const real G1_3_4 = map.det*c0_3*c1_4 + map.det*c1_3*c0_4;
  const real G1_3_5 = map.det*c0_3*c1_5 + map.det*c1_3*c0_5;
  const real G1_4_1 = map.det*c0_4*c1_1 + map.det*c1_4*c0_1;
  const real G1_4_3 = map.det*c0_4*c1_3 + map.det*c1_4*c0_3;
  const real G1_4_4 = map.det*c0_4*c1_4 + map.det*c1_4*c0_4;
  const real G1_4_5 = map.det*c0_4*c1_5 + map.det*c1_4*c0_5;
  const real G1_5_2 = map.det*c0_5*c1_2 + map.det*c1_5*c0_2;
  const real G1_5_3 = map.det*c0_5*c1_3 + map.det*c1_5*c0_3;
  const real G1_5_4 = map.det*c0_5*c1_4 + map.det*c1_5*c0_4;
  const real G1_5_5 = map.det*c0_5*c1_5 + map.det*c1_5*c0_5;
  const real G1_6_6 = map.det*c0_6*c1_6 + map.det*c1_6*c0_6;
  const real G1_6_7 = map.det*c0_6*c1_7 + map.det*c1_6*c0_7;
  const real G1_6_8 = map.det*c0_6*c1_8 + map.det*c1_6*c0_8;
  const real G1_6_9 = map.det*c0_6*c1_9 + map.det*c1_6*c0_9;
  const real G1_7_6 = map.det*c0_7*c1_6 + map.det*c1_7*c0_6;
  const real G1_7_7 = map.det*c0_7*c1_7 + map.det*c1_7*c0_7;
  const real G1_7_8 = map.det*c0_7*c1_8 + map.det*c1_7*c0_8;
  const real G1_7_10 = map.det*c0_7*c1_10 + map.det*c1_7*c0_10;
  const real G1_8_6 = map.det*c0_8*c1_6 + map.det*c1_8*c0_6;
  const real G1_8_7 = map.det*c0_8*c1_7 + map.det*c1_8*c0_7;
  const real G1_8_8 = map.det*c0_8*c1_8 + map.det*c1_8*c0_8;
  const real G1_8_11 = map.det*c0_8*c1_11 + map.det*c1_8*c0_11;
  const real G1_9_6 = map.det*c0_9*c1_6 + map.det*c1_9*c0_6;
  const real G1_9_9 = map.det*c0_9*c1_9 + map.det*c1_9*c0_9;
  const real G1_9_10 = map.det*c0_9*c1_10 + map.det*c1_9*c0_10;
  const real G1_9_11 = map.det*c0_9*c1_11 + map.det*c1_9*c0_11;
  const real G1_10_7 = map.det*c0_10*c1_7 + map.det*c1_10*c0_7;
  const real G1_10_9 = map.det*c0_10*c1_9 + map.det*c1_10*c0_9;
  const real G1_10_10 = map.det*c0_10*c1_10 + map.det*c1_10*c0_10;
  const real G1_10_11 = map.det*c0_10*c1_11 + map.det*c1_10*c0_11;
  const real G1_11_8 = map.det*c0_11*c1_8 + map.det*c1_11*c0_8;
  const real G1_11_9 = map.det*c0_11*c1_9 + map.det*c1_11*c0_9;
  const real G1_11_10 = map.det*c0_11*c1_10 + map.det*c1_11*c0_10;
  const real G1_11_11 = map.det*c0_11*c1_11 + map.det*c1_11*c0_11;

  // Compute element tensor
  block[0] = 1.666666666666665e-02*G0_0_0 - 2.777777777777774e-03*G0_0_1 - 2.777777777777775e-03*G0_0_2 - 1.111111111111110e-02*G0_0_3 - 2.777777777777774e-03*G0_1_0 + 1.666666666666665e-02*G0_1_1 - 2.777777777777776e-03*G0_1_2 - 1.111111111111111e-02*G0_1_4 - 2.777777777777775e-03*G0_2_0 - 2.777777777777776e-03*G0_2_1 + 1.666666666666666e-02*G0_2_2 - 1.111111111111111e-02*G0_2_5 - 1.111111111111110e-02*G0_3_0 + 8.888888888888882e-02*G0_3_3 + 4.444444444444443e-02*G0_3_4 + 4.444444444444443e-02*G0_3_5 - 1.111111111111111e-02*G0_4_1 + 4.444444444444443e-02*G0_4_3 + 8.888888888888884e-02*G0_4_4 + 4.444444444444442e-02*G0_4_5 - 1.111111111111111e-02*G0_5_2 + 4.444444444444443e-02*G0_5_3 + 4.444444444444443e-02*G0_5_4 + 8.888888888888882e-02*G0_5_5 + 1.666666666666665e-02*G0_6_6 - 2.777777777777774e-03*G0_6_7 - 2.777777777777774e-03*G0_6_8 - 1.111111111111109e-02*G0_6_9 - 2.777777777777774e-03*G0_7_6 + 1.666666666666665e-02*G0_7_7 - 2.777777777777775e-03*G0_7_8 - 1.111111111111111e-02*G0_7_10 - 2.777777777777774e-03*G0_8_6 - 2.777777777777775e-03*G0_8_7 + 1.666666666666666e-02*G0_8_8 - 1.111111111111111e-02*G0_8_11 - 1.111111111111109e-02*G0_9_6 + 8.888888888888882e-02*G0_9_9 + 4.444444444444443e-02*G0_9_10 + 4.444444444444443e-02*G0_9_11 - 1.111111111111111e-02*G0_10_7 + 4.444444444444443e-02*G0_10_9 + 8.888888888888884e-02*G0_10_10 + 4.444444444444442e-02*G0_10_11 - 1.111111111111111e-02*G0_11_8 + 4.444444444444443e-02*G0_11_9 + 4.444444444444443e-02*G0_11_10 + 8.888888888888882e-02*G0_11_11 - 1.666666666666665e-02*G1_0_0 + 2.777777777777774e-03*G1_0_1 + 2.777777777777775e-03*G1_0_2 + 1.111111111111110e-02*G1_0_3 + 2.777777777777774e-03*G1_1_0 - 1.666666666666665e-02*G1_1_1 + 2.777777777777776e-03*G1_1_2 + 1.111111111111111e-02*G1_1_4 + 2.777777777777775e-03*G1_2_0 + 2.777777777777776e-03*G1_2_1 - 1.666666666666666e-02*G1_2_2 + 1.111111111111111e-02*G1_2_5 + 1.111111111111110e-02*G1_3_0 - 8.888888888888882e-02*G1_3_3 - 4.444444444444443e-02*G1_3_4 - 4.444444444444443e-02*G1_3_5 + 1.111111111111111e-02*G1_4_1 - 4.444444444444443e-02*G1_4_3 - 8.888888888888884e-02*G1_4_4 - 4.444444444444442e-02*G1_4_5 + 1.111111111111111e-02*G1_5_2 - 4.444444444444443e-02*G1_5_3 - 4.444444444444443e-02*G1_5_4 - 8.888888888888882e-02*G1_5_5 - 1.666666666666665e-02*G1_6_6 + 2.777777777777774e-03*G1_6_7 + 2.777777777777774e-03*G1_6_8 + 1.111111111111109e-02*G1_6_9 + 2.777777777777774e-03*G1_7_6 - 1.666666666666665e-02*G1_7_7 + 2.777777777777775e-03*G1_7_8 + 1.111111111111111e-02*G1_7_10 + 2.777777777777774e-03*G1_8_6 + 2.777777777777775e-03*G1_8_7 - 1.666666666666666e-02*G1_8_8 + 1.111111111111111e-02*G1_8_11 + 1.111111111111109e-02*G1_9_6 - 8.888888888888882e-02*G1_9_9 - 4.444444444444443e-02*G1_9_10 - 4.444444444444443e-02*G1_9_11 + 1.111111111111111e-02*G1_10_7 - 4.444444444444443e-02*G1_10_9 - 8.888888888888884e-02*G1_10_10 - 4.444444444444442e-02*G1_10_11 + 1.111111111111111e-02*G1_11_8 - 4.444444444444443e-02*G1_11_9 - 4.444444444444443e-02*G1_11_10 - 8.888888888888882e-02*G1_11_11;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}

// No contribution from interior boundaries
bool LinearForm::interior_boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
