// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.3-dev.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __BIDOMAIN_H
#define __BIDOMAIN_H

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

namespace dolfin { namespace Bidomain {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  class FunctionElement_0;

  class FunctionElement_1;

  class FunctionElement_2;

  class FunctionElement_3;

  class FunctionElement_4;

  class FunctionElement_5;

  class FunctionElement_6;

  class FunctionElement_7;

  class FunctionElement_8;

  BilinearForm(Function& w0, Function& w1, Function& w2, Function& w3, Function& w4, Function& w5, Function& w6, Function& w7, Function& w8);
  

  void eval(real block[], const AffineMap& map) const;

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
    nodes[0] = cell.vertexID(0);
    nodes[1] = cell.vertexID(1);
    nodes[2] = cell.vertexID(2);
    nodes[3] = cell.vertexID(3);
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
    nodes[0] = cell.vertexID(0);
    nodes[1] = cell.vertexID(1);
    nodes[2] = cell.vertexID(2);
    nodes[3] = cell.vertexID(3);
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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_1 : public dolfin::FiniteElement
{
public:

  FunctionElement_1() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_2 : public dolfin::FiniteElement
{
public:

  FunctionElement_2() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_3 : public dolfin::FiniteElement
{
public:

  FunctionElement_3() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_4 : public dolfin::FiniteElement
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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_5 : public dolfin::FiniteElement
{
public:

  FunctionElement_5() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_5()
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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_6 : public dolfin::FiniteElement
{
public:

  FunctionElement_6() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_6()
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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_7 : public dolfin::FiniteElement
{
public:

  FunctionElement_7() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_7()
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
    nodes[0] = cell.id();
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

class BilinearForm::FunctionElement_8 : public dolfin::FiniteElement
{
public:

  FunctionElement_8() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~FunctionElement_8()
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
    nodes[0] = cell.id();
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

BilinearForm::BilinearForm(Function& w0, Function& w1, Function& w2, Function& w3, Function& w4, Function& w5, Function& w6, Function& w7, Function& w8) : dolfin::BilinearForm(9)
{
  // Create finite element for test space
  _test = new TestElement();

  // Create finite element for trial space
  _trial = new TrialElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
  initFunction(1, w1, new FunctionElement_1());
  initFunction(2, w2, new FunctionElement_2());
  initFunction(3, w3, new FunctionElement_3());
  initFunction(4, w4, new FunctionElement_4());
  initFunction(5, w5, new FunctionElement_5());
  initFunction(6, w6, new FunctionElement_6());
  initFunction(7, w7, new FunctionElement_7());
  initFunction(8, w8, new FunctionElement_8());
}

void BilinearForm::eval(real block[], const AffineMap& map) const
{
  // Compute coefficients
  const real c0_0 = c[0][0];
  const real c1_0 = c[1][0];
  const real c2_0 = c[2][0];
  const real c3_0 = c[3][0];
  const real c4_0 = c[4][0];
  const real c5_0 = c[5][0];
  const real c6_0 = c[6][0];
  const real c7_0 = c[7][0];
  const real c8_0 = c[8][0];

  // Compute geometry tensors
  const real G0_0_0_0 = map.det*c0_0*map.g00*map.g00 + map.det*c1_0*map.g00*map.g01 + map.det*c2_0*map.g00*map.g02 + map.det*c3_0*map.g01*map.g00 + map.det*c4_0*map.g01*map.g01 + map.det*c5_0*map.g01*map.g02 + map.det*c6_0*map.g02*map.g00 + map.det*c7_0*map.g02*map.g01 + map.det*c8_0*map.g02*map.g02;
  const real G0_0_0_1 = map.det*c0_0*map.g00*map.g10 + map.det*c1_0*map.g00*map.g11 + map.det*c2_0*map.g00*map.g12 + map.det*c3_0*map.g01*map.g10 + map.det*c4_0*map.g01*map.g11 + map.det*c5_0*map.g01*map.g12 + map.det*c6_0*map.g02*map.g10 + map.det*c7_0*map.g02*map.g11 + map.det*c8_0*map.g02*map.g12;
  const real G0_0_0_2 = map.det*c0_0*map.g00*map.g20 + map.det*c1_0*map.g00*map.g21 + map.det*c2_0*map.g00*map.g22 + map.det*c3_0*map.g01*map.g20 + map.det*c4_0*map.g01*map.g21 + map.det*c5_0*map.g01*map.g22 + map.det*c6_0*map.g02*map.g20 + map.det*c7_0*map.g02*map.g21 + map.det*c8_0*map.g02*map.g22;
  const real G0_1_0_0 = map.det*c0_0*map.g10*map.g00 + map.det*c1_0*map.g10*map.g01 + map.det*c2_0*map.g10*map.g02 + map.det*c3_0*map.g11*map.g00 + map.det*c4_0*map.g11*map.g01 + map.det*c5_0*map.g11*map.g02 + map.det*c6_0*map.g12*map.g00 + map.det*c7_0*map.g12*map.g01 + map.det*c8_0*map.g12*map.g02;
  const real G0_1_0_1 = map.det*c0_0*map.g10*map.g10 + map.det*c1_0*map.g10*map.g11 + map.det*c2_0*map.g10*map.g12 + map.det*c3_0*map.g11*map.g10 + map.det*c4_0*map.g11*map.g11 + map.det*c5_0*map.g11*map.g12 + map.det*c6_0*map.g12*map.g10 + map.det*c7_0*map.g12*map.g11 + map.det*c8_0*map.g12*map.g12;
  const real G0_1_0_2 = map.det*c0_0*map.g10*map.g20 + map.det*c1_0*map.g10*map.g21 + map.det*c2_0*map.g10*map.g22 + map.det*c3_0*map.g11*map.g20 + map.det*c4_0*map.g11*map.g21 + map.det*c5_0*map.g11*map.g22 + map.det*c6_0*map.g12*map.g20 + map.det*c7_0*map.g12*map.g21 + map.det*c8_0*map.g12*map.g22;
  const real G0_2_0_0 = map.det*c0_0*map.g20*map.g00 + map.det*c1_0*map.g20*map.g01 + map.det*c2_0*map.g20*map.g02 + map.det*c3_0*map.g21*map.g00 + map.det*c4_0*map.g21*map.g01 + map.det*c5_0*map.g21*map.g02 + map.det*c6_0*map.g22*map.g00 + map.det*c7_0*map.g22*map.g01 + map.det*c8_0*map.g22*map.g02;
  const real G0_2_0_1 = map.det*c0_0*map.g20*map.g10 + map.det*c1_0*map.g20*map.g11 + map.det*c2_0*map.g20*map.g12 + map.det*c3_0*map.g21*map.g10 + map.det*c4_0*map.g21*map.g11 + map.det*c5_0*map.g21*map.g12 + map.det*c6_0*map.g22*map.g10 + map.det*c7_0*map.g22*map.g11 + map.det*c8_0*map.g22*map.g12;
  const real G0_2_0_2 = map.det*c0_0*map.g20*map.g20 + map.det*c1_0*map.g20*map.g21 + map.det*c2_0*map.g20*map.g22 + map.det*c3_0*map.g21*map.g20 + map.det*c4_0*map.g21*map.g21 + map.det*c5_0*map.g21*map.g22 + map.det*c6_0*map.g22*map.g20 + map.det*c7_0*map.g22*map.g21 + map.det*c8_0*map.g22*map.g22;

  // Compute element tensor
  block[0] = 1.666666666666664e-01*G0_0_0_0 + 1.666666666666664e-01*G0_0_0_1 + 1.666666666666664e-01*G0_0_0_2 + 1.666666666666664e-01*G0_1_0_0 + 1.666666666666665e-01*G0_1_0_1 + 1.666666666666665e-01*G0_1_0_2 + 1.666666666666664e-01*G0_2_0_0 + 1.666666666666665e-01*G0_2_0_1 + 1.666666666666665e-01*G0_2_0_2;
  block[1] = -1.666666666666664e-01*G0_0_0_0 - 1.666666666666664e-01*G0_1_0_0 - 1.666666666666664e-01*G0_2_0_0;
  block[2] = -1.666666666666664e-01*G0_0_0_1 - 1.666666666666665e-01*G0_1_0_1 - 1.666666666666665e-01*G0_2_0_1;
  block[3] = -1.666666666666664e-01*G0_0_0_2 - 1.666666666666665e-01*G0_1_0_2 - 1.666666666666665e-01*G0_2_0_2;
  block[4] = -1.666666666666664e-01*G0_0_0_0 - 1.666666666666664e-01*G0_0_0_1 - 1.666666666666664e-01*G0_0_0_2;
  block[5] = 1.666666666666664e-01*G0_0_0_0;
  block[6] = 1.666666666666664e-01*G0_0_0_1;
  block[7] = 1.666666666666664e-01*G0_0_0_2;
  block[8] = -1.666666666666664e-01*G0_1_0_0 - 1.666666666666665e-01*G0_1_0_1 - 1.666666666666665e-01*G0_1_0_2;
  block[9] = 1.666666666666664e-01*G0_1_0_0;
  block[10] = 1.666666666666665e-01*G0_1_0_1;
  block[11] = 1.666666666666665e-01*G0_1_0_2;
  block[12] = -1.666666666666664e-01*G0_2_0_0 - 1.666666666666665e-01*G0_2_0_1 - 1.666666666666665e-01*G0_2_0_2;
  block[13] = 1.666666666666664e-01*G0_2_0_0;
  block[14] = 1.666666666666665e-01*G0_2_0_1;
  block[15] = 1.666666666666665e-01*G0_2_0_2;
}

// No contribution from the boundary
void BilinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}   
} }

#endif
