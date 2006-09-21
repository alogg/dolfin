// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.3-dev.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __POISSON_H
#define __POISSON_H

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

namespace dolfin { namespace Poisson {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  BilinearForm();
  

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
    nodes[0] = cell.vertexID(0);
    nodes[1] = cell.vertexID(1);
    nodes[2] = cell.vertexID(2);
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
    nodes[0] = cell.vertexID(0);
    nodes[1] = cell.vertexID(1);
    nodes[2] = cell.vertexID(2);
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

BilinearForm::BilinearForm() : dolfin::BilinearForm(0)
{
  // Create finite element for test space
  _test = new TestElement();

  // Create finite element for trial space
  _trial = new TrialElement();
}

void BilinearForm::eval(real block[], const AffineMap& map) const
{
  // Compute geometry tensors
  const real G0_0_0 = map.det*map.g00*map.g00 + map.det*map.g01*map.g01;
  const real G0_0_1 = map.det*map.g00*map.g10 + map.det*map.g01*map.g11;
  const real G0_1_0 = map.det*map.g10*map.g00 + map.det*map.g11*map.g01;
  const real G0_1_1 = map.det*map.g10*map.g10 + map.det*map.g11*map.g11;

  // Compute element tensor
  block[0] = 4.999999999999998e-01*G0_0_0 + 4.999999999999997e-01*G0_0_1 + 4.999999999999997e-01*G0_1_0 + 4.999999999999996e-01*G0_1_1;
  block[1] = -4.999999999999998e-01*G0_0_0 - 4.999999999999997e-01*G0_1_0;
  block[2] = -4.999999999999997e-01*G0_0_1 - 4.999999999999996e-01*G0_1_1;
  block[3] = -4.999999999999998e-01*G0_0_0 - 4.999999999999997e-01*G0_0_1;
  block[4] = 4.999999999999998e-01*G0_0_0;
  block[5] = 4.999999999999997e-01*G0_0_1;
  block[6] = -4.999999999999997e-01*G0_1_0 - 4.999999999999996e-01*G0_1_1;
  block[7] = 4.999999999999997e-01*G0_1_0;
  block[8] = 4.999999999999996e-01*G0_1_1;
}

// No contribution from the boundary
void BilinearForm::eval(real block[], const AffineMap& map, unsigned int facet) const {}   
} }

#endif
