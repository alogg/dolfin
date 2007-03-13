// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
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
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

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
    nodes[0] = 3*cell.index() + 0;
    nodes[1] = 3*cell.index() + 1;
    nodes[2] = 3*cell.index() + 2;
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
    FiniteElementSpec s("Discontinuous Lagrange", "triangle", 1);
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
    nodes[0] = 3*cell.index() + 0;
    nodes[1] = 3*cell.index() + 1;
    nodes[2] = 3*cell.index() + 2;
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
    FiniteElementSpec s("Discontinuous Lagrange", "triangle", 1);
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

// No contribution from the interior
bool BilinearForm::interior_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map, real det) const {}

// No contribution from the boundary
bool BilinearForm::boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}


// Contribution from interior boundaries
bool BilinearForm::interior_boundary_contribution() const { return true; }

void BilinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const
{
  // Compute geometry tensors
  const real G0_ = det;

  // Compute interior facet tensor
  switch ( facet0 )
  {
  case 0:
    switch ( facet1 )
    {
    case 0:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 0.000000000000000e+00;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 3.333333333333331e-01*G0_;
        block[11] = 1.666666666666666e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 1.666666666666665e-01*G0_;
        block[17] = 3.333333333333330e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 0.000000000000000e+00;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 1.666666666666666e-01*G0_;
        block[11] = 3.333333333333330e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 3.333333333333330e-01*G0_;
        block[17] = 1.666666666666665e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    case 1:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 0.000000000000000e+00;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 1.666666666666666e-01*G0_;
        block[10] = 0.000000000000000e+00;
        block[11] = 3.333333333333330e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 3.333333333333330e-01*G0_;
        block[16] = 0.000000000000000e+00;
        block[17] = 1.666666666666665e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 0.000000000000000e+00;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 3.333333333333330e-01*G0_;
        block[10] = 0.000000000000000e+00;
        block[11] = 1.666666666666666e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 1.666666666666665e-01*G0_;
        block[16] = 0.000000000000000e+00;
        block[17] = 3.333333333333330e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    case 2:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 0.000000000000000e+00;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 3.333333333333330e-01*G0_;
        block[10] = 1.666666666666666e-01*G0_;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 1.666666666666665e-01*G0_;
        block[16] = 3.333333333333330e-01*G0_;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 0.000000000000000e+00;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 1.666666666666665e-01*G0_;
        block[10] = 3.333333333333331e-01*G0_;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 3.333333333333330e-01*G0_;
        block[16] = 1.666666666666665e-01*G0_;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    }
    break;
  case 1:
    switch ( facet1 )
    {
    case 0:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 1.666666666666665e-01*G0_;
        block[5] = 3.333333333333330e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 0.000000000000000e+00;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 3.333333333333330e-01*G0_;
        block[17] = 1.666666666666665e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 3.333333333333330e-01*G0_;
        block[5] = 1.666666666666665e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 0.000000000000000e+00;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 1.666666666666665e-01*G0_;
        block[17] = 3.333333333333330e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    case 1:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 3.333333333333330e-01*G0_;
        block[4] = 0.000000000000000e+00;
        block[5] = 1.666666666666665e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 0.000000000000000e+00;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 1.666666666666665e-01*G0_;
        block[16] = 0.000000000000000e+00;
        block[17] = 3.333333333333330e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 1.666666666666665e-01*G0_;
        block[4] = 0.000000000000000e+00;
        block[5] = 3.333333333333330e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 0.000000000000000e+00;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 3.333333333333330e-01*G0_;
        block[16] = 0.000000000000000e+00;
        block[17] = 1.666666666666665e-01*G0_;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    case 2:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 1.666666666666665e-01*G0_;
        block[4] = 3.333333333333330e-01*G0_;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 0.000000000000000e+00;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 3.333333333333330e-01*G0_;
        block[16] = 1.666666666666665e-01*G0_;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 3.333333333333330e-01*G0_;
        block[4] = 1.666666666666665e-01*G0_;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 0.000000000000000e+00;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 1.666666666666665e-01*G0_;
        block[16] = 3.333333333333330e-01*G0_;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    }
    break;
  case 2:
    switch ( facet1 )
    {
    case 0:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 3.333333333333330e-01*G0_;
        block[5] = 1.666666666666665e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 1.666666666666666e-01*G0_;
        block[11] = 3.333333333333330e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 0.000000000000000e+00;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 0.000000000000000e+00;
        block[4] = 1.666666666666665e-01*G0_;
        block[5] = 3.333333333333330e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 0.000000000000000e+00;
        block[10] = 3.333333333333331e-01*G0_;
        block[11] = 1.666666666666665e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 0.000000000000000e+00;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    case 1:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 1.666666666666665e-01*G0_;
        block[4] = 0.000000000000000e+00;
        block[5] = 3.333333333333330e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 3.333333333333330e-01*G0_;
        block[10] = 0.000000000000000e+00;
        block[11] = 1.666666666666665e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 0.000000000000000e+00;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 3.333333333333330e-01*G0_;
        block[4] = 0.000000000000000e+00;
        block[5] = 1.666666666666665e-01*G0_;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 1.666666666666665e-01*G0_;
        block[10] = 0.000000000000000e+00;
        block[11] = 3.333333333333330e-01*G0_;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 0.000000000000000e+00;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    case 2:
      switch ( alignment )
      {
      case 0:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 3.333333333333330e-01*G0_;
        block[4] = 1.666666666666665e-01*G0_;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 1.666666666666665e-01*G0_;
        block[10] = 3.333333333333331e-01*G0_;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 0.000000000000000e+00;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      case 1:
        block[0] = 0.000000000000000e+00;
        block[1] = 0.000000000000000e+00;
        block[2] = 0.000000000000000e+00;
        block[3] = 1.666666666666665e-01*G0_;
        block[4] = 3.333333333333330e-01*G0_;
        block[5] = 0.000000000000000e+00;
        block[6] = 0.000000000000000e+00;
        block[7] = 0.000000000000000e+00;
        block[8] = 0.000000000000000e+00;
        block[9] = 3.333333333333330e-01*G0_;
        block[10] = 1.666666666666665e-01*G0_;
        block[11] = 0.000000000000000e+00;
        block[12] = 0.000000000000000e+00;
        block[13] = 0.000000000000000e+00;
        block[14] = 0.000000000000000e+00;
        block[15] = 0.000000000000000e+00;
        block[16] = 0.000000000000000e+00;
        block[17] = 0.000000000000000e+00;
        block[18] = 0.000000000000000e+00;
        block[19] = 0.000000000000000e+00;
        block[20] = 0.000000000000000e+00;
        block[21] = 0.000000000000000e+00;
        block[22] = 0.000000000000000e+00;
        block[23] = 0.000000000000000e+00;
        block[24] = 0.000000000000000e+00;
        block[25] = 0.000000000000000e+00;
        block[26] = 0.000000000000000e+00;
        block[27] = 0.000000000000000e+00;
        block[28] = 0.000000000000000e+00;
        block[29] = 0.000000000000000e+00;
        block[30] = 0.000000000000000e+00;
        block[31] = 0.000000000000000e+00;
        block[32] = 0.000000000000000e+00;
        block[33] = 0.000000000000000e+00;
        block[34] = 0.000000000000000e+00;
        block[35] = 0.000000000000000e+00;
        break;
      }
      break;
    }
    break;
  }
}

} }

#endif
