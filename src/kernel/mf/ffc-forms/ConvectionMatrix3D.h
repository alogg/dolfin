// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.3-dev.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __CONVECTIONMATRIX3D_H
#define __CONVECTIONMATRIX3D_H

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

namespace dolfin { namespace ConvectionMatrix3D {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  BilinearForm(const real& c0, const real& c1, const real& c2);
  

  void eval(real block[], const AffineMap& map) const;

  void eval(real block[], const AffineMap& map, unsigned int facet) const;

private:

  const real& c0;  const real& c1;  const real& c2;

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

BilinearForm::BilinearForm(const real& c0, const real& c1, const real& c2) : dolfin::BilinearForm(0), c0(c0), c1(c1), c2(c2)
{
  // Create finite element for test space
  _test = new TestElement();

  // Create finite element for trial space
  _trial = new TrialElement();
}

void BilinearForm::eval(real block[], const AffineMap& map) const
{
  // Compute geometry tensors
  const real G0_0 = map.det*c0*map.g00 + map.det*c1*map.g01 + map.det*c2*map.g02;
  const real G0_1 = map.det*c0*map.g10 + map.det*c1*map.g11 + map.det*c2*map.g12;
  const real G0_2 = map.det*c0*map.g20 + map.det*c1*map.g21 + map.det*c2*map.g22;

  // Compute element tensor
  block[0] = -4.166666666666661e-02*G0_0 - 4.166666666666662e-02*G0_1 - 4.166666666666662e-02*G0_2;
  block[1] = 4.166666666666661e-02*G0_0;
  block[2] = 4.166666666666662e-02*G0_1;
  block[3] = 4.166666666666662e-02*G0_2;
  block[4] = -4.166666666666661e-02*G0_0 - 4.166666666666662e-02*G0_1 - 4.166666666666662e-02*G0_2;
  block[5] = 4.166666666666661e-02*G0_0;
  block[6] = 4.166666666666662e-02*G0_1;
  block[7] = 4.166666666666662e-02*G0_2;
  block[8] = -4.166666666666661e-02*G0_0 - 4.166666666666662e-02*G0_1 - 4.166666666666662e-02*G0_2;
  block[9] = 4.166666666666661e-02*G0_0;
  block[10] = 4.166666666666662e-02*G0_1;
  block[11] = 4.166666666666662e-02*G0_2;
  block[12] = -4.166666666666661e-02*G0_0 - 4.166666666666662e-02*G0_1 - 4.166666666666662e-02*G0_2;
  block[13] = 4.166666666666661e-02*G0_0;
  block[14] = 4.166666666666662e-02*G0_1;
  block[15] = 4.166666666666662e-02*G0_2;
}

// No contribution from the boundary
void BilinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}   
} }

#endif
