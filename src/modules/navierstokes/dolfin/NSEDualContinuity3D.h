// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.4.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __NSEDUALCONTINUITY3D_H
#define __NSEDUALCONTINUITY3D_H

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

namespace dolfin { namespace NSEDualContinuity3D {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  class FunctionElement_0;

  BilinearForm(Function& w0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

};

class BilinearForm::TestElement : public dolfin::FiniteElement
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
    return 4;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    nodes[3] = cell.entities(0)[3];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 1);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class BilinearForm::TrialElement : public dolfin::FiniteElement
{
public:

  TrialElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~TrialElement()
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
    return 4;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    nodes[3] = cell.entities(0)[3];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 1);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class BilinearForm::FunctionElement_0 : public dolfin::FiniteElement
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
    return 1;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
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
    FiniteElementSpec s("Discontinuous Lagrange", "tetrahedron", 0);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

BilinearForm::BilinearForm(Function& w0) : dolfin::BilinearForm(1)
{
  // Create finite element for test space
  _test = new TestElement();

  // Create finite element for trial space
  _trial = new TrialElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
}

// Contribution from the interior
bool BilinearForm::interior_contribution() const { return true; }

void BilinearForm::eval(real block[], const AffineMap& map) const
{
  // Compute coefficients
  const real c0_0 = c[0][0];

  // Compute geometry tensors
  const real G0_0_0_0 = map.det*c0_0*map.g00*map.g00 + map.det*c0_0*map.g01*map.g01 + map.det*c0_0*map.g02*map.g02;
  const real G0_0_0_1 = map.det*c0_0*map.g00*map.g10 + map.det*c0_0*map.g01*map.g11 + map.det*c0_0*map.g02*map.g12;
  const real G0_0_0_2 = map.det*c0_0*map.g00*map.g20 + map.det*c0_0*map.g01*map.g21 + map.det*c0_0*map.g02*map.g22;
  const real G0_0_1_0 = map.det*c0_0*map.g10*map.g00 + map.det*c0_0*map.g11*map.g01 + map.det*c0_0*map.g12*map.g02;
  const real G0_0_1_1 = map.det*c0_0*map.g10*map.g10 + map.det*c0_0*map.g11*map.g11 + map.det*c0_0*map.g12*map.g12;
  const real G0_0_1_2 = map.det*c0_0*map.g10*map.g20 + map.det*c0_0*map.g11*map.g21 + map.det*c0_0*map.g12*map.g22;
  const real G0_0_2_0 = map.det*c0_0*map.g20*map.g00 + map.det*c0_0*map.g21*map.g01 + map.det*c0_0*map.g22*map.g02;
  const real G0_0_2_1 = map.det*c0_0*map.g20*map.g10 + map.det*c0_0*map.g21*map.g11 + map.det*c0_0*map.g22*map.g12;
  const real G0_0_2_2 = map.det*c0_0*map.g20*map.g20 + map.det*c0_0*map.g21*map.g21 + map.det*c0_0*map.g22*map.g22;

  // Compute element tensor
  block[0] = 1.666666666666665e-01*G0_0_0_0 + 1.666666666666665e-01*G0_0_0_1 + 1.666666666666664e-01*G0_0_0_2 + 1.666666666666665e-01*G0_0_1_0 + 1.666666666666665e-01*G0_0_1_1 + 1.666666666666664e-01*G0_0_1_2 + 1.666666666666664e-01*G0_0_2_0 + 1.666666666666664e-01*G0_0_2_1 + 1.666666666666664e-01*G0_0_2_2;
  block[1] = -1.666666666666665e-01*G0_0_0_0 - 1.666666666666665e-01*G0_0_1_0 - 1.666666666666664e-01*G0_0_2_0;
  block[2] = -1.666666666666665e-01*G0_0_0_1 - 1.666666666666665e-01*G0_0_1_1 - 1.666666666666664e-01*G0_0_2_1;
  block[3] = -1.666666666666665e-01*G0_0_0_2 - 1.666666666666665e-01*G0_0_1_2 - 1.666666666666664e-01*G0_0_2_2;
  block[4] = -1.666666666666665e-01*G0_0_0_0 - 1.666666666666665e-01*G0_0_0_1 - 1.666666666666664e-01*G0_0_0_2;
  block[5] = 1.666666666666665e-01*G0_0_0_0;
  block[6] = 1.666666666666665e-01*G0_0_0_1;
  block[7] = 1.666666666666665e-01*G0_0_0_2;
  block[8] = -1.666666666666665e-01*G0_0_1_0 - 1.666666666666665e-01*G0_0_1_1 - 1.666666666666664e-01*G0_0_1_2;
  block[9] = 1.666666666666665e-01*G0_0_1_0;
  block[10] = 1.666666666666665e-01*G0_0_1_1;
  block[11] = 1.666666666666665e-01*G0_0_1_2;
  block[12] = -1.666666666666665e-01*G0_0_2_0 - 1.666666666666665e-01*G0_0_2_1 - 1.666666666666664e-01*G0_0_2_2;
  block[13] = 1.666666666666665e-01*G0_0_2_0;
  block[14] = 1.666666666666665e-01*G0_0_2_1;
  block[15] = 1.666666666666665e-01*G0_0_2_2;
}

// No contribution from the boundary
bool BilinearForm::boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:

  class TestElement;

  class FunctionElement_0;

  LinearForm(Function& w0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

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
    return 4;
  }

  inline unsigned int shapedim() const
  {
    return 3;
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
    nodes[3] = cell.entities(0)[3];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 1);
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

LinearForm::LinearForm(Function& w0) : dolfin::LinearForm(1)
{
  // Create finite element for test space
  _test = new TestElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
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

  // Compute geometry tensors
  const real G0_0_0_0 = map.det*c0_0*map.g00;
  const real G0_0_1_0 = map.det*c0_0*map.g10;
  const real G0_0_2_0 = map.det*c0_0*map.g20;
  const real G0_1_0_0 = map.det*c0_1*map.g00;
  const real G0_2_1_0 = map.det*c0_2*map.g10;
  const real G0_3_2_0 = map.det*c0_3*map.g20;
  const real G0_4_0_1 = map.det*c0_4*map.g01;
  const real G0_4_1_1 = map.det*c0_4*map.g11;
  const real G0_4_2_1 = map.det*c0_4*map.g21;
  const real G0_5_0_1 = map.det*c0_5*map.g01;
  const real G0_6_1_1 = map.det*c0_6*map.g11;
  const real G0_7_2_1 = map.det*c0_7*map.g21;
  const real G0_8_0_2 = map.det*c0_8*map.g02;
  const real G0_8_1_2 = map.det*c0_8*map.g12;
  const real G0_8_2_2 = map.det*c0_8*map.g22;
  const real G0_9_0_2 = map.det*c0_9*map.g02;
  const real G0_10_1_2 = map.det*c0_10*map.g12;
  const real G0_11_2_2 = map.det*c0_11*map.g22;

  // Compute element tensor
  block[0] = 4.166666666666662e-02*G0_0_0_0 + 4.166666666666662e-02*G0_0_1_0 + 4.166666666666662e-02*G0_0_2_0 - 4.166666666666662e-02*G0_1_0_0 - 4.166666666666662e-02*G0_2_1_0 - 4.166666666666662e-02*G0_3_2_0 + 4.166666666666662e-02*G0_4_0_1 + 4.166666666666662e-02*G0_4_1_1 + 4.166666666666662e-02*G0_4_2_1 - 4.166666666666662e-02*G0_5_0_1 - 4.166666666666662e-02*G0_6_1_1 - 4.166666666666662e-02*G0_7_2_1 + 4.166666666666662e-02*G0_8_0_2 + 4.166666666666662e-02*G0_8_1_2 + 4.166666666666662e-02*G0_8_2_2 - 4.166666666666662e-02*G0_9_0_2 - 4.166666666666662e-02*G0_10_1_2 - 4.166666666666662e-02*G0_11_2_2;
  block[1] = 4.166666666666662e-02*G0_0_0_0 + 4.166666666666662e-02*G0_0_1_0 + 4.166666666666661e-02*G0_0_2_0 - 4.166666666666662e-02*G0_1_0_0 - 4.166666666666662e-02*G0_2_1_0 - 4.166666666666662e-02*G0_3_2_0 + 4.166666666666662e-02*G0_4_0_1 + 4.166666666666662e-02*G0_4_1_1 + 4.166666666666661e-02*G0_4_2_1 - 4.166666666666662e-02*G0_5_0_1 - 4.166666666666662e-02*G0_6_1_1 - 4.166666666666662e-02*G0_7_2_1 + 4.166666666666662e-02*G0_8_0_2 + 4.166666666666662e-02*G0_8_1_2 + 4.166666666666661e-02*G0_8_2_2 - 4.166666666666662e-02*G0_9_0_2 - 4.166666666666662e-02*G0_10_1_2 - 4.166666666666662e-02*G0_11_2_2;
  block[2] = 4.166666666666662e-02*G0_0_0_0 + 4.166666666666662e-02*G0_0_1_0 + 4.166666666666661e-02*G0_0_2_0 - 4.166666666666662e-02*G0_1_0_0 - 4.166666666666662e-02*G0_2_1_0 - 4.166666666666662e-02*G0_3_2_0 + 4.166666666666662e-02*G0_4_0_1 + 4.166666666666662e-02*G0_4_1_1 + 4.166666666666661e-02*G0_4_2_1 - 4.166666666666662e-02*G0_5_0_1 - 4.166666666666662e-02*G0_6_1_1 - 4.166666666666662e-02*G0_7_2_1 + 4.166666666666662e-02*G0_8_0_2 + 4.166666666666662e-02*G0_8_1_2 + 4.166666666666661e-02*G0_8_2_2 - 4.166666666666662e-02*G0_9_0_2 - 4.166666666666662e-02*G0_10_1_2 - 4.166666666666662e-02*G0_11_2_2;
  block[3] = 4.166666666666662e-02*G0_0_0_0 + 4.166666666666662e-02*G0_0_1_0 + 4.166666666666661e-02*G0_0_2_0 - 4.166666666666662e-02*G0_1_0_0 - 4.166666666666662e-02*G0_2_1_0 - 4.166666666666662e-02*G0_3_2_0 + 4.166666666666662e-02*G0_4_0_1 + 4.166666666666662e-02*G0_4_1_1 + 4.166666666666661e-02*G0_4_2_1 - 4.166666666666662e-02*G0_5_0_1 - 4.166666666666662e-02*G0_6_1_1 - 4.166666666666662e-02*G0_7_2_1 + 4.166666666666662e-02*G0_8_0_2 + 4.166666666666662e-02*G0_8_1_2 + 4.166666666666661e-02*G0_8_2_2 - 4.166666666666662e-02*G0_9_0_2 - 4.166666666666662e-02*G0_10_1_2 - 4.166666666666662e-02*G0_11_2_2;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}

} }

#endif
