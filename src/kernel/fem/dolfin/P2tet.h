// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.3.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __P2TET_H
#define __P2TET_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>

namespace dolfin
{
  
  class P2tet : public dolfin::FiniteElement
  {
  public:
  
    P2tet() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~P2tet()
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
      return 10;
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
      int offset = mesh.numVertices();
      nodes[4] = offset + cell.edgeID(0);
      nodes[5] = offset + cell.edgeID(1);
      nodes[6] = offset + cell.edgeID(2);
      nodes[7] = offset + cell.edgeID(3);
      nodes[8] = offset + cell.edgeID(4);
      nodes[9] = offset + cell.edgeID(5);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      points[4] = map(5.000000000000000e-01, 5.000000000000000e-01, 0.000000000000000e+00);
      points[5] = map(0.000000000000000e+00, 5.000000000000000e-01, 0.000000000000000e+00);
      points[6] = map(5.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
      points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 5.000000000000000e-01);
      points[8] = map(5.000000000000000e-01, 0.000000000000000e+00, 5.000000000000000e-01);
      points[9] = map(0.000000000000000e+00, 5.000000000000000e-01, 5.000000000000000e-01);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
      components[6] = 0;
      components[7] = 0;
      components[8] = 0;
      components[9] = 0;
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
      FiniteElementSpec s("Lagrange", "tetrahedron", 2);
      return s;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
  };
  
}

#endif
