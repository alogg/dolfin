// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __STIFFNESSMATRIX3D_H
#define __STIFFNESSMATRIX3D_H

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

namespace dolfin { namespace StiffnessMatrix3D {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  BilinearForm(const real& c0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

private:

  const real& c0;

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
    error("Element is scalar.");
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
    error("Element is scalar.");
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

BilinearForm::BilinearForm(const real& c0) : dolfin::BilinearForm(0), c0(c0)
{
  // Create finite element for test space
  _test = new TestElement();

  // Create finite element for trial space
  _trial = new TrialElement();
}

// Contribution from the interior
bool BilinearForm::interior_contribution() const { return true; }

void BilinearForm::eval(real block[], const AffineMap& map, real det) const
{
  // Compute geometry tensors
  const real G0_0_0 = det*c0*map.g00*map.g00 + det*c0*map.g01*map.g01 + det*c0*map.g02*map.g02;
  const real G0_0_1 = det*c0*map.g00*map.g10 + det*c0*map.g01*map.g11 + det*c0*map.g02*map.g12;
  const real G0_0_2 = det*c0*map.g00*map.g20 + det*c0*map.g01*map.g21 + det*c0*map.g02*map.g22;
  const real G0_1_0 = det*c0*map.g10*map.g00 + det*c0*map.g11*map.g01 + det*c0*map.g12*map.g02;
  const real G0_1_1 = det*c0*map.g10*map.g10 + det*c0*map.g11*map.g11 + det*c0*map.g12*map.g12;
  const real G0_1_2 = det*c0*map.g10*map.g20 + det*c0*map.g11*map.g21 + det*c0*map.g12*map.g22;
  const real G0_2_0 = det*c0*map.g20*map.g00 + det*c0*map.g21*map.g01 + det*c0*map.g22*map.g02;
  const real G0_2_1 = det*c0*map.g20*map.g10 + det*c0*map.g21*map.g11 + det*c0*map.g22*map.g12;
  const real G0_2_2 = det*c0*map.g20*map.g20 + det*c0*map.g21*map.g21 + det*c0*map.g22*map.g22;

  // Compute element tensor
  block[0] = 1.666666666666664e-01*G0_0_0 + 1.666666666666664e-01*G0_0_1 + 1.666666666666664e-01*G0_0_2 + 1.666666666666664e-01*G0_1_0 + 1.666666666666665e-01*G0_1_1 + 1.666666666666665e-01*G0_1_2 + 1.666666666666664e-01*G0_2_0 + 1.666666666666665e-01*G0_2_1 + 1.666666666666665e-01*G0_2_2;
  block[1] = -1.666666666666664e-01*G0_0_0 - 1.666666666666664e-01*G0_1_0 - 1.666666666666664e-01*G0_2_0;
  block[2] = -1.666666666666664e-01*G0_0_1 - 1.666666666666665e-01*G0_1_1 - 1.666666666666665e-01*G0_2_1;
  block[3] = -1.666666666666664e-01*G0_0_2 - 1.666666666666665e-01*G0_1_2 - 1.666666666666665e-01*G0_2_2;
  block[4] = -1.666666666666664e-01*G0_0_0 - 1.666666666666664e-01*G0_0_1 - 1.666666666666664e-01*G0_0_2;
  block[5] = 1.666666666666664e-01*G0_0_0;
  block[6] = 1.666666666666664e-01*G0_0_1;
  block[7] = 1.666666666666664e-01*G0_0_2;
  block[8] = -1.666666666666664e-01*G0_1_0 - 1.666666666666665e-01*G0_1_1 - 1.666666666666665e-01*G0_1_2;
  block[9] = 1.666666666666664e-01*G0_1_0;
  block[10] = 1.666666666666665e-01*G0_1_1;
  block[11] = 1.666666666666665e-01*G0_1_2;
  block[12] = -1.666666666666664e-01*G0_2_0 - 1.666666666666665e-01*G0_2_1 - 1.666666666666665e-01*G0_2_2;
  block[13] = 1.666666666666664e-01*G0_2_0;
  block[14] = 1.666666666666665e-01*G0_2_1;
  block[15] = 1.666666666666665e-01*G0_2_2;
}

// No contribution from the boundary
bool BilinearForm::boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool BilinearForm::interior_boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
