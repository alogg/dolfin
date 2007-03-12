// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __CNSRESMOMENTUM_H
#define __CNSRESMOMENTUM_H

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

namespace dolfin { namespace CNSResMomentum {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:

  class TestElement;

  class FunctionElement_0;

  class FunctionElement_1;

  class FunctionElement_2;

  class FunctionElement_3;

  class FunctionElement_4;

  LinearForm(Function& w0, Function& w1, Function& w2, Function& w3, Function& w4, const real& c0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

private:

  const real& c0;

};

class LinearForm::TestElement : public dolfin::FiniteElement
{
public:

  TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 2;

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
    return 2;
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
    nodes[0] = cell.index();
    int offset = mesh.topology().size(2);
    nodes[1] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(3.333333333333334e-01, 3.333333333333334e-01);
    points[1] = map(3.333333333333334e-01, 3.333333333333334e-01);
    components[0] = 0;
    components[1] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.numCells();
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
    FiniteElementSpec s("Discontinuous vector Lagrange", "triangle", 0, 2);
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
    return 3;
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
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
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
    FiniteElementSpec s("Lagrange", "triangle", 1);
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
    return 6;
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
    nodes[3] = offset + cell.entities(0)[0];
    nodes[4] = offset + cell.entities(0)[1];
    nodes[5] = offset + cell.entities(0)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[4] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 1;
    components[4] = 1;
    components[5] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 1, 2);
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
    tensordims[0] = 2;

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
    return 6;
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
    nodes[3] = offset + cell.entities(0)[0];
    nodes[4] = offset + cell.entities(0)[1];
    nodes[5] = offset + cell.entities(0)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[4] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 1;
    components[4] = 1;
    components[5] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 1, 2);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_3 : public dolfin::FiniteElement
{
public:

  FunctionElement_3() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 2;

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_3()
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
    nodes[3] = offset + cell.entities(0)[0];
    nodes[4] = offset + cell.entities(0)[1];
    nodes[5] = offset + cell.entities(0)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[4] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 1;
    components[4] = 1;
    components[5] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 1, 2);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_4 : public dolfin::FiniteElement
{
public:

  FunctionElement_4() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_4()
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

LinearForm::LinearForm(Function& w0, Function& w1, Function& w2, Function& w3, Function& w4, const real& c0) : dolfin::LinearForm(5), c0(c0)
{
  // Create finite element for test space
  _test = new TestElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
  initFunction(1, w1, new FunctionElement_1());
  initFunction(2, w2, new FunctionElement_2());
  initFunction(3, w3, new FunctionElement_3());
  initFunction(4, w4, new FunctionElement_4());
}

// Contribution from the interior
bool LinearForm::interior_contribution() const { return true; }

void LinearForm::eval(real block[], const AffineMap& map, real det) const
{
  // Compute coefficients
  const real c0_0 = c[0][0];
  const real c0_1 = c[0][1];
  const real c0_2 = c[0][2];
  const real c1_0 = c[1][0];
  const real c1_1 = c[1][1];
  const real c1_2 = c[1][2];
  const real c1_3 = c[1][3];
  const real c1_4 = c[1][4];
  const real c1_5 = c[1][5];
  const real c2_0 = c[2][0];
  const real c2_1 = c[2][1];
  const real c2_2 = c[2][2];
  const real c2_3 = c[2][3];
  const real c2_4 = c[2][4];
  const real c2_5 = c[2][5];
  const real c3_0 = c[3][0];
  const real c3_1 = c[3][1];
  const real c3_2 = c[3][2];
  const real c3_3 = c[3][3];
  const real c3_4 = c[3][4];
  const real c3_5 = c[3][5];
  const real c4_0 = c[4][0];

  // Compute geometry tensors
  const real G0_0_0 = det*(1.0/c0)*c2_0*c4_0;
  const real G0_1_0 = det*(1.0/c0)*c2_1*c4_0;
  const real G0_2_0 = det*(1.0/c0)*c2_2*c4_0;
  const real G0_3_0 = det*(1.0/c0)*c2_3*c4_0;
  const real G0_4_0 = det*(1.0/c0)*c2_4*c4_0;
  const real G0_5_0 = det*(1.0/c0)*c2_5*c4_0;
  const real G1_0_0 = det*(1.0/c0)*c3_0*c4_0;
  const real G1_1_0 = det*(1.0/c0)*c3_1*c4_0;
  const real G1_2_0 = det*(1.0/c0)*c3_2*c4_0;
  const real G1_3_0 = det*(1.0/c0)*c3_3*c4_0;
  const real G1_4_0 = det*(1.0/c0)*c3_4*c4_0;
  const real G1_5_0 = det*(1.0/c0)*c3_5*c4_0;
  const real G2_0_0_0_0 = det*c1_0*c2_0*c4_0*map.g00;
  const real G2_0_0_1_0 = det*c1_0*c2_0*c4_0*map.g10;
  const real G2_0_1_0_0 = det*c1_0*c2_1*c4_0*map.g00;
  const real G2_0_1_1_0 = det*c1_0*c2_1*c4_0*map.g10;
  const real G2_0_2_0_0 = det*c1_0*c2_2*c4_0*map.g00;
  const real G2_0_2_1_0 = det*c1_0*c2_2*c4_0*map.g10;
  const real G2_1_0_0_0 = det*c1_1*c2_0*c4_0*map.g00;
  const real G2_1_1_0_0 = det*c1_1*c2_1*c4_0*map.g00;
  const real G2_1_2_0_0 = det*c1_1*c2_2*c4_0*map.g00;
  const real G2_2_0_1_0 = det*c1_2*c2_0*c4_0*map.g10;
  const real G2_2_1_1_0 = det*c1_2*c2_1*c4_0*map.g10;
  const real G2_2_2_1_0 = det*c1_2*c2_2*c4_0*map.g10;
  const real G3_0_0_0_0 = det*c1_0*c2_0*c4_0*map.g00;
  const real G3_0_0_1_0 = det*c1_0*c2_0*c4_0*map.g10;
  const real G3_0_1_0_0 = det*c1_0*c2_1*c4_0*map.g00;
  const real G3_0_2_1_0 = det*c1_0*c2_2*c4_0*map.g10;
  const real G3_1_0_0_0 = det*c1_1*c2_0*c4_0*map.g00;
  const real G3_1_0_1_0 = det*c1_1*c2_0*c4_0*map.g10;
  const real G3_1_1_0_0 = det*c1_1*c2_1*c4_0*map.g00;
  const real G3_1_2_1_0 = det*c1_1*c2_2*c4_0*map.g10;
  const real G3_2_0_0_0 = det*c1_2*c2_0*c4_0*map.g00;
  const real G3_2_0_1_0 = det*c1_2*c2_0*c4_0*map.g10;
  const real G3_2_1_0_0 = det*c1_2*c2_1*c4_0*map.g00;
  const real G3_2_2_1_0 = det*c1_2*c2_2*c4_0*map.g10;
  const real G4_3_0_0_0 = det*c1_3*c2_0*c4_0*map.g01;
  const real G4_3_0_1_0 = det*c1_3*c2_0*c4_0*map.g11;
  const real G4_3_1_0_0 = det*c1_3*c2_1*c4_0*map.g01;
  const real G4_3_1_1_0 = det*c1_3*c2_1*c4_0*map.g11;
  const real G4_3_2_0_0 = det*c1_3*c2_2*c4_0*map.g01;
  const real G4_3_2_1_0 = det*c1_3*c2_2*c4_0*map.g11;
  const real G4_4_0_0_0 = det*c1_4*c2_0*c4_0*map.g01;
  const real G4_4_1_0_0 = det*c1_4*c2_1*c4_0*map.g01;
  const real G4_4_2_0_0 = det*c1_4*c2_2*c4_0*map.g01;
  const real G4_5_0_1_0 = det*c1_5*c2_0*c4_0*map.g11;
  const real G4_5_1_1_0 = det*c1_5*c2_1*c4_0*map.g11;
  const real G4_5_2_1_0 = det*c1_5*c2_2*c4_0*map.g11;
  const real G5_3_0_0_0 = det*c1_3*c2_0*c4_0*map.g01;
  const real G5_3_0_1_0 = det*c1_3*c2_0*c4_0*map.g11;
  const real G5_3_1_0_0 = det*c1_3*c2_1*c4_0*map.g01;
  const real G5_3_2_1_0 = det*c1_3*c2_2*c4_0*map.g11;
  const real G5_4_0_0_0 = det*c1_4*c2_0*c4_0*map.g01;
  const real G5_4_0_1_0 = det*c1_4*c2_0*c4_0*map.g11;
  const real G5_4_1_0_0 = det*c1_4*c2_1*c4_0*map.g01;
  const real G5_4_2_1_0 = det*c1_4*c2_2*c4_0*map.g11;
  const real G5_5_0_0_0 = det*c1_5*c2_0*c4_0*map.g01;
  const real G5_5_0_1_0 = det*c1_5*c2_0*c4_0*map.g11;
  const real G5_5_1_0_0 = det*c1_5*c2_1*c4_0*map.g01;
  const real G5_5_2_1_0 = det*c1_5*c2_2*c4_0*map.g11;
  const real G6_0_3_0_0 = det*c1_0*c2_3*c4_0*map.g00;
  const real G6_0_3_1_0 = det*c1_0*c2_3*c4_0*map.g10;
  const real G6_0_4_0_0 = det*c1_0*c2_4*c4_0*map.g00;
  const real G6_0_4_1_0 = det*c1_0*c2_4*c4_0*map.g10;
  const real G6_0_5_0_0 = det*c1_0*c2_5*c4_0*map.g00;
  const real G6_0_5_1_0 = det*c1_0*c2_5*c4_0*map.g10;
  const real G6_1_3_0_0 = det*c1_1*c2_3*c4_0*map.g00;
  const real G6_1_4_0_0 = det*c1_1*c2_4*c4_0*map.g00;
  const real G6_1_5_0_0 = det*c1_1*c2_5*c4_0*map.g00;
  const real G6_2_3_1_0 = det*c1_2*c2_3*c4_0*map.g10;
  const real G6_2_4_1_0 = det*c1_2*c2_4*c4_0*map.g10;
  const real G6_2_5_1_0 = det*c1_2*c2_5*c4_0*map.g10;
  const real G7_0_3_0_0 = det*c1_0*c2_3*c4_0*map.g00;
  const real G7_0_3_1_0 = det*c1_0*c2_3*c4_0*map.g10;
  const real G7_0_4_0_0 = det*c1_0*c2_4*c4_0*map.g00;
  const real G7_0_5_1_0 = det*c1_0*c2_5*c4_0*map.g10;
  const real G7_1_3_0_0 = det*c1_1*c2_3*c4_0*map.g00;
  const real G7_1_3_1_0 = det*c1_1*c2_3*c4_0*map.g10;
  const real G7_1_4_0_0 = det*c1_1*c2_4*c4_0*map.g00;
  const real G7_1_5_1_0 = det*c1_1*c2_5*c4_0*map.g10;
  const real G7_2_3_0_0 = det*c1_2*c2_3*c4_0*map.g00;
  const real G7_2_3_1_0 = det*c1_2*c2_3*c4_0*map.g10;
  const real G7_2_4_0_0 = det*c1_2*c2_4*c4_0*map.g00;
  const real G7_2_5_1_0 = det*c1_2*c2_5*c4_0*map.g10;
  const real G8_3_3_0_0 = det*c1_3*c2_3*c4_0*map.g01;
  const real G8_3_3_1_0 = det*c1_3*c2_3*c4_0*map.g11;
  const real G8_3_4_0_0 = det*c1_3*c2_4*c4_0*map.g01;
  const real G8_3_4_1_0 = det*c1_3*c2_4*c4_0*map.g11;
  const real G8_3_5_0_0 = det*c1_3*c2_5*c4_0*map.g01;
  const real G8_3_5_1_0 = det*c1_3*c2_5*c4_0*map.g11;
  const real G8_4_3_0_0 = det*c1_4*c2_3*c4_0*map.g01;
  const real G8_4_4_0_0 = det*c1_4*c2_4*c4_0*map.g01;
  const real G8_4_5_0_0 = det*c1_4*c2_5*c4_0*map.g01;
  const real G8_5_3_1_0 = det*c1_5*c2_3*c4_0*map.g11;
  const real G8_5_4_1_0 = det*c1_5*c2_4*c4_0*map.g11;
  const real G8_5_5_1_0 = det*c1_5*c2_5*c4_0*map.g11;
  const real G9_3_3_0_0 = det*c1_3*c2_3*c4_0*map.g01;
  const real G9_3_3_1_0 = det*c1_3*c2_3*c4_0*map.g11;
  const real G9_3_4_0_0 = det*c1_3*c2_4*c4_0*map.g01;
  const real G9_3_5_1_0 = det*c1_3*c2_5*c4_0*map.g11;
  const real G9_4_3_0_0 = det*c1_4*c2_3*c4_0*map.g01;
  const real G9_4_3_1_0 = det*c1_4*c2_3*c4_0*map.g11;
  const real G9_4_4_0_0 = det*c1_4*c2_4*c4_0*map.g01;
  const real G9_4_5_1_0 = det*c1_4*c2_5*c4_0*map.g11;
  const real G9_5_3_0_0 = det*c1_5*c2_3*c4_0*map.g01;
  const real G9_5_3_1_0 = det*c1_5*c2_3*c4_0*map.g11;
  const real G9_5_4_0_0 = det*c1_5*c2_4*c4_0*map.g01;
  const real G9_5_5_1_0 = det*c1_5*c2_5*c4_0*map.g11;
  const real G10_0_0_0 = det*c0_0*c4_0*map.g00;
  const real G10_0_1_0 = det*c0_0*c4_0*map.g10;
  const real G10_1_0_0 = det*c0_1*c4_0*map.g00;
  const real G10_2_1_0 = det*c0_2*c4_0*map.g10;
  const real G11_0_0_0 = det*c0_0*c4_0*map.g01;
  const real G11_0_1_0 = det*c0_0*c4_0*map.g11;
  const real G11_1_0_0 = det*c0_1*c4_0*map.g01;
  const real G11_2_1_0 = det*c0_2*c4_0*map.g11;

  // Compute element tensor
  block[0] = 1.666666666666665e-01*G0_0_0 + 1.666666666666665e-01*G0_1_0 + 1.666666666666665e-01*G0_2_0 - 1.666666666666665e-01*G1_0_0 - 1.666666666666665e-01*G1_1_0 - 1.666666666666665e-01*G1_2_0 - 1.666666666666665e-01*G2_0_0_0_0 - 1.666666666666665e-01*G2_0_0_1_0 - 1.666666666666666e-01*G2_0_1_0_0 - 1.666666666666665e-01*G2_0_1_1_0 - 1.666666666666665e-01*G2_0_2_0_0 - 1.666666666666665e-01*G2_0_2_1_0 + 1.666666666666665e-01*G2_1_0_0_0 + 1.666666666666666e-01*G2_1_1_0_0 + 1.666666666666665e-01*G2_1_2_0_0 + 1.666666666666665e-01*G2_2_0_1_0 + 1.666666666666665e-01*G2_2_1_1_0 + 1.666666666666665e-01*G2_2_2_1_0 - 1.666666666666665e-01*G3_0_0_0_0 - 1.666666666666665e-01*G3_0_0_1_0 + 1.666666666666665e-01*G3_0_1_0_0 + 1.666666666666665e-01*G3_0_2_1_0 - 1.666666666666666e-01*G3_1_0_0_0 - 1.666666666666665e-01*G3_1_0_1_0 + 1.666666666666666e-01*G3_1_1_0_0 + 1.666666666666665e-01*G3_1_2_1_0 - 1.666666666666665e-01*G3_2_0_0_0 - 1.666666666666665e-01*G3_2_0_1_0 + 1.666666666666665e-01*G3_2_1_0_0 + 1.666666666666665e-01*G3_2_2_1_0 - 1.666666666666665e-01*G4_3_0_0_0 - 1.666666666666665e-01*G4_3_0_1_0 - 1.666666666666666e-01*G4_3_1_0_0 - 1.666666666666665e-01*G4_3_1_1_0 - 1.666666666666665e-01*G4_3_2_0_0 - 1.666666666666665e-01*G4_3_2_1_0 + 1.666666666666665e-01*G4_4_0_0_0 + 1.666666666666666e-01*G4_4_1_0_0 + 1.666666666666665e-01*G4_4_2_0_0 + 1.666666666666665e-01*G4_5_0_1_0 + 1.666666666666665e-01*G4_5_1_1_0 + 1.666666666666665e-01*G4_5_2_1_0 - 1.666666666666665e-01*G5_3_0_0_0 - 1.666666666666665e-01*G5_3_0_1_0 + 1.666666666666665e-01*G5_3_1_0_0 + 1.666666666666665e-01*G5_3_2_1_0 - 1.666666666666666e-01*G5_4_0_0_0 - 1.666666666666665e-01*G5_4_0_1_0 + 1.666666666666666e-01*G5_4_1_0_0 + 1.666666666666665e-01*G5_4_2_1_0 - 1.666666666666665e-01*G5_5_0_0_0 - 1.666666666666665e-01*G5_5_0_1_0 + 1.666666666666665e-01*G5_5_1_0_0 + 1.666666666666665e-01*G5_5_2_1_0 - 4.999999999999997e-01*G10_0_0_0 - 4.999999999999996e-01*G10_0_1_0 + 4.999999999999997e-01*G10_1_0_0 + 4.999999999999996e-01*G10_2_1_0;
  block[1] = 1.666666666666665e-01*G0_3_0 + 1.666666666666665e-01*G0_4_0 + 1.666666666666665e-01*G0_5_0 - 1.666666666666665e-01*G1_3_0 - 1.666666666666665e-01*G1_4_0 - 1.666666666666665e-01*G1_5_0 - 1.666666666666665e-01*G6_0_3_0_0 - 1.666666666666665e-01*G6_0_3_1_0 - 1.666666666666666e-01*G6_0_4_0_0 - 1.666666666666665e-01*G6_0_4_1_0 - 1.666666666666665e-01*G6_0_5_0_0 - 1.666666666666665e-01*G6_0_5_1_0 + 1.666666666666665e-01*G6_1_3_0_0 + 1.666666666666666e-01*G6_1_4_0_0 + 1.666666666666665e-01*G6_1_5_0_0 + 1.666666666666665e-01*G6_2_3_1_0 + 1.666666666666665e-01*G6_2_4_1_0 + 1.666666666666665e-01*G6_2_5_1_0 - 1.666666666666665e-01*G7_0_3_0_0 - 1.666666666666665e-01*G7_0_3_1_0 + 1.666666666666665e-01*G7_0_4_0_0 + 1.666666666666665e-01*G7_0_5_1_0 - 1.666666666666666e-01*G7_1_3_0_0 - 1.666666666666665e-01*G7_1_3_1_0 + 1.666666666666666e-01*G7_1_4_0_0 + 1.666666666666665e-01*G7_1_5_1_0 - 1.666666666666665e-01*G7_2_3_0_0 - 1.666666666666665e-01*G7_2_3_1_0 + 1.666666666666665e-01*G7_2_4_0_0 + 1.666666666666665e-01*G7_2_5_1_0 - 1.666666666666665e-01*G8_3_3_0_0 - 1.666666666666665e-01*G8_3_3_1_0 - 1.666666666666666e-01*G8_3_4_0_0 - 1.666666666666665e-01*G8_3_4_1_0 - 1.666666666666665e-01*G8_3_5_0_0 - 1.666666666666665e-01*G8_3_5_1_0 + 1.666666666666665e-01*G8_4_3_0_0 + 1.666666666666666e-01*G8_4_4_0_0 + 1.666666666666665e-01*G8_4_5_0_0 + 1.666666666666665e-01*G8_5_3_1_0 + 1.666666666666665e-01*G8_5_4_1_0 + 1.666666666666665e-01*G8_5_5_1_0 - 1.666666666666665e-01*G9_3_3_0_0 - 1.666666666666665e-01*G9_3_3_1_0 + 1.666666666666665e-01*G9_3_4_0_0 + 1.666666666666665e-01*G9_3_5_1_0 - 1.666666666666666e-01*G9_4_3_0_0 - 1.666666666666665e-01*G9_4_3_1_0 + 1.666666666666666e-01*G9_4_4_0_0 + 1.666666666666665e-01*G9_4_5_1_0 - 1.666666666666665e-01*G9_5_3_0_0 - 1.666666666666665e-01*G9_5_3_1_0 + 1.666666666666665e-01*G9_5_4_0_0 + 1.666666666666665e-01*G9_5_5_1_0 - 4.999999999999997e-01*G11_0_0_0 - 4.999999999999996e-01*G11_0_1_0 + 4.999999999999997e-01*G11_1_0_0 + 4.999999999999996e-01*G11_2_1_0;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool LinearForm::interior_boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
