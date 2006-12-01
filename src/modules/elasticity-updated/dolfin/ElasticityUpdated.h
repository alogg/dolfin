// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.4.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __ELASTICITYUPDATED_H
#define __ELASTICITYUPDATED_H

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

namespace dolfin { namespace ElasticityUpdated {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:

  class TestElement;

  class FunctionElement_0;

  class FunctionElement_1;

  class FunctionElement_2;

  LinearForm(Function& w0, Function& w1, Function& w2, const real& c0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

private:

  const real& c0;

};

class LinearForm::TestElement : public dolfin::FiniteElement
{
public:

  TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 3;

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
    return 12;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(0)[0];
    nodes[5] = offset + cell.entities(0)[1];
    nodes[6] = offset + cell.entities(0)[2];
    nodes[7] = offset + cell.entities(0)[3];
    offset = offset + mesh.topology().size(0);
    nodes[8] = offset + cell.entities(0)[0];
    nodes[9] = offset + cell.entities(0)[1];
    nodes[10] = offset + cell.entities(0)[2];
    nodes[11] = offset + cell.entities(0)[3];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[6] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[8] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[9] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[10] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[11] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 1;
    components[5] = 1;
    components[6] = 1;
    components[7] = 1;
    components[8] = 2;
    components[9] = 2;
    components[10] = 2;
    components[11] = 2;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(0);
    vertex_nodes[2] = offset + vertex;
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
    FiniteElementSpec s("Vector Lagrange", "tetrahedron", 1, 3);
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
    tensordims[0] = 3;

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
    return 3;
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
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(0)[0];
    nodes[5] = offset + cell.entities(0)[1];
    nodes[6] = offset + cell.entities(0)[2];
    nodes[7] = offset + cell.entities(0)[3];
    offset = offset + mesh.topology().size(0);
    nodes[8] = offset + cell.entities(0)[0];
    nodes[9] = offset + cell.entities(0)[1];
    nodes[10] = offset + cell.entities(0)[2];
    nodes[11] = offset + cell.entities(0)[3];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[6] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[8] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[9] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[10] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[11] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 1;
    components[5] = 1;
    components[6] = 1;
    components[7] = 1;
    components[8] = 2;
    components[9] = 2;
    components[10] = 2;
    components[11] = 2;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(0);
    vertex_nodes[2] = offset + vertex;
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
    FiniteElementSpec s("Vector Lagrange", "tetrahedron", 1, 3);
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
    tensordims[0] = 9;

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
    return 9;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    nodes[0] = cell.index();
    int offset = mesh.topology().size(3);
    nodes[1] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[2] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[3] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[4] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[5] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[6] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[7] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[8] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[1] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[2] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[3] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[4] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[5] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[6] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[7] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[8] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    components[0] = 0;
    components[1] = 1;
    components[2] = 2;
    components[3] = 3;
    components[4] = 4;
    components[5] = 5;
    components[6] = 6;
    components[7] = 7;
    components[8] = 8;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(3);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[2] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[3] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[4] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[5] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[6] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[7] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[8] = offset + vertex;
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
    FiniteElementSpec s("Discontinuous vector Lagrange", "tetrahedron", 0, 9);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_2 : public dolfin::FiniteElement
{
public:

  FunctionElement_2() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 9;

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_2()
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
    return 9;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    nodes[0] = cell.index();
    int offset = mesh.topology().size(3);
    nodes[1] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[2] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[3] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[4] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[5] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[6] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[7] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[8] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[1] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[2] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[3] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[4] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[5] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[6] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[7] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[8] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    components[0] = 0;
    components[1] = 1;
    components[2] = 2;
    components[3] = 3;
    components[4] = 4;
    components[5] = 5;
    components[6] = 6;
    components[7] = 7;
    components[8] = 8;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(3);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[2] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[3] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[4] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[5] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[6] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[7] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[8] = offset + vertex;
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
    FiniteElementSpec s("Discontinuous vector Lagrange", "tetrahedron", 0, 9);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

LinearForm::LinearForm(Function& w0, Function& w1, Function& w2, const real& c0) : dolfin::LinearForm(3), c0(c0)
{
  // Create finite element for test space
  _test = new TestElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
  initFunction(1, w1, new FunctionElement_1());
  initFunction(2, w2, new FunctionElement_2());
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
  const real c2_0 = c[2][0];
  const real c2_1 = c[2][1];
  const real c2_2 = c[2][2];
  const real c2_3 = c[2][3];
  const real c2_4 = c[2][4];
  const real c2_5 = c[2][5];
  const real c2_6 = c[2][6];
  const real c2_7 = c[2][7];
  const real c2_8 = c[2][8];

  // Compute geometry tensors
  const real G0_0 = map.det*c0_0;
  const real G0_1 = map.det*c0_1;
  const real G0_2 = map.det*c0_2;
  const real G0_3 = map.det*c0_3;
  const real G0_4 = map.det*c0_4;
  const real G0_5 = map.det*c0_5;
  const real G0_6 = map.det*c0_6;
  const real G0_7 = map.det*c0_7;
  const real G0_8 = map.det*c0_8;
  const real G0_9 = map.det*c0_9;
  const real G0_10 = map.det*c0_10;
  const real G0_11 = map.det*c0_11;
  const real G1_0_0 = map.det*c1_0*map.g00 + map.det*c0*c2_0*map.g00;
  const real G1_0_1 = map.det*c1_0*map.g10 + map.det*c0*c2_0*map.g10;
  const real G1_0_2 = map.det*c1_0*map.g20 + map.det*c0*c2_0*map.g20;
  const real G2_3_0 = map.det*c1_3*map.g01 + map.det*c0*c2_3*map.g01;
  const real G2_3_1 = map.det*c1_3*map.g11 + map.det*c0*c2_3*map.g11;
  const real G2_3_2 = map.det*c1_3*map.g21 + map.det*c0*c2_3*map.g21;
  const real G3_6_0 = map.det*c1_6*map.g02 + map.det*c0*c2_6*map.g02;
  const real G3_6_1 = map.det*c1_6*map.g12 + map.det*c0*c2_6*map.g12;
  const real G3_6_2 = map.det*c1_6*map.g22 + map.det*c0*c2_6*map.g22;
  const real G4_1_0 = map.det*c1_1*map.g00 + map.det*c0*c2_1*map.g00;
  const real G4_1_1 = map.det*c1_1*map.g10 + map.det*c0*c2_1*map.g10;
  const real G4_1_2 = map.det*c1_1*map.g20 + map.det*c0*c2_1*map.g20;
  const real G5_4_0 = map.det*c1_4*map.g01 + map.det*c0*c2_4*map.g01;
  const real G5_4_1 = map.det*c1_4*map.g11 + map.det*c0*c2_4*map.g11;
  const real G5_4_2 = map.det*c1_4*map.g21 + map.det*c0*c2_4*map.g21;
  const real G6_7_0 = map.det*c1_7*map.g02 + map.det*c0*c2_7*map.g02;
  const real G6_7_1 = map.det*c1_7*map.g12 + map.det*c0*c2_7*map.g12;
  const real G6_7_2 = map.det*c1_7*map.g22 + map.det*c0*c2_7*map.g22;
  const real G7_2_0 = map.det*c1_2*map.g00 + map.det*c0*c2_2*map.g00;
  const real G7_2_1 = map.det*c1_2*map.g10 + map.det*c0*c2_2*map.g10;
  const real G7_2_2 = map.det*c1_2*map.g20 + map.det*c0*c2_2*map.g20;
  const real G8_5_0 = map.det*c1_5*map.g01 + map.det*c0*c2_5*map.g01;
  const real G8_5_1 = map.det*c1_5*map.g11 + map.det*c0*c2_5*map.g11;
  const real G8_5_2 = map.det*c1_5*map.g21 + map.det*c0*c2_5*map.g21;
  const real G9_8_0 = map.det*c1_8*map.g02 + map.det*c0*c2_8*map.g02;
  const real G9_8_1 = map.det*c1_8*map.g12 + map.det*c0*c2_8*map.g12;
  const real G9_8_2 = map.det*c1_8*map.g22 + map.det*c0*c2_8*map.g22;

  // Compute element tensor
  block[0] = 1.666666666666662e-02*G0_0 + 8.333333333333309e-03*G0_1 + 8.333333333333309e-03*G0_2 + 8.333333333333311e-03*G0_3 + 1.666666666666664e-01*G1_0_0 + 1.666666666666665e-01*G1_0_1 + 1.666666666666664e-01*G1_0_2 + 1.666666666666664e-01*G2_3_0 + 1.666666666666665e-01*G2_3_1 + 1.666666666666664e-01*G2_3_2 + 1.666666666666664e-01*G3_6_0 + 1.666666666666665e-01*G3_6_1 + 1.666666666666664e-01*G3_6_2;
  block[1] = 8.333333333333307e-03*G0_0 + 1.666666666666661e-02*G0_1 + 8.333333333333309e-03*G0_2 + 8.333333333333309e-03*G0_3 - 1.666666666666664e-01*G1_0_0 - 1.666666666666664e-01*G2_3_0 - 1.666666666666664e-01*G3_6_0;
  block[2] = 8.333333333333309e-03*G0_0 + 8.333333333333309e-03*G0_1 + 1.666666666666662e-02*G0_2 + 8.333333333333311e-03*G0_3 - 1.666666666666665e-01*G1_0_1 - 1.666666666666665e-01*G2_3_1 - 1.666666666666665e-01*G3_6_1;
  block[3] = 8.333333333333311e-03*G0_0 + 8.333333333333307e-03*G0_1 + 8.333333333333311e-03*G0_2 + 1.666666666666662e-02*G0_3 - 1.666666666666665e-01*G1_0_2 - 1.666666666666665e-01*G2_3_2 - 1.666666666666665e-01*G3_6_2;
  block[4] = 1.666666666666662e-02*G0_4 + 8.333333333333307e-03*G0_5 + 8.333333333333309e-03*G0_6 + 8.333333333333309e-03*G0_7 + 1.666666666666664e-01*G4_1_0 + 1.666666666666665e-01*G4_1_1 + 1.666666666666665e-01*G4_1_2 + 1.666666666666664e-01*G5_4_0 + 1.666666666666665e-01*G5_4_1 + 1.666666666666665e-01*G5_4_2 + 1.666666666666664e-01*G6_7_0 + 1.666666666666665e-01*G6_7_1 + 1.666666666666665e-01*G6_7_2;
  block[5] = 8.333333333333307e-03*G0_4 + 1.666666666666662e-02*G0_5 + 8.333333333333312e-03*G0_6 + 8.333333333333311e-03*G0_7 - 1.666666666666664e-01*G4_1_0 - 1.666666666666664e-01*G5_4_0 - 1.666666666666664e-01*G6_7_0;
  block[6] = 8.333333333333311e-03*G0_4 + 8.333333333333312e-03*G0_5 + 1.666666666666662e-02*G0_6 + 8.333333333333314e-03*G0_7 - 1.666666666666665e-01*G4_1_1 - 1.666666666666665e-01*G5_4_1 - 1.666666666666665e-01*G6_7_1;
  block[7] = 8.333333333333309e-03*G0_4 + 8.333333333333311e-03*G0_5 + 8.333333333333314e-03*G0_6 + 1.666666666666662e-02*G0_7 - 1.666666666666665e-01*G4_1_2 - 1.666666666666665e-01*G5_4_2 - 1.666666666666665e-01*G6_7_2;
  block[8] = 1.666666666666662e-02*G0_8 + 8.333333333333309e-03*G0_9 + 8.333333333333307e-03*G0_10 + 8.333333333333309e-03*G0_11 + 1.666666666666664e-01*G7_2_0 + 1.666666666666665e-01*G7_2_1 + 1.666666666666665e-01*G7_2_2 + 1.666666666666664e-01*G8_5_0 + 1.666666666666665e-01*G8_5_1 + 1.666666666666665e-01*G8_5_2 + 1.666666666666664e-01*G9_8_0 + 1.666666666666665e-01*G9_8_1 + 1.666666666666665e-01*G9_8_2;
  block[9] = 8.333333333333309e-03*G0_8 + 1.666666666666662e-02*G0_9 + 8.333333333333311e-03*G0_10 + 8.333333333333312e-03*G0_11 - 1.666666666666664e-01*G7_2_0 - 1.666666666666664e-01*G8_5_0 - 1.666666666666664e-01*G9_8_0;
  block[10] = 8.333333333333307e-03*G0_8 + 8.333333333333311e-03*G0_9 + 1.666666666666662e-02*G0_10 + 8.333333333333311e-03*G0_11 - 1.666666666666665e-01*G7_2_1 - 1.666666666666665e-01*G8_5_1 - 1.666666666666665e-01*G9_8_1;
  block[11] = 8.333333333333309e-03*G0_8 + 8.333333333333311e-03*G0_9 + 8.333333333333311e-03*G0_10 + 1.666666666666662e-02*G0_11 - 1.666666666666665e-01*G7_2_2 - 1.666666666666665e-01*G8_5_2 - 1.666666666666665e-01*G9_8_2;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}

} }

#endif
