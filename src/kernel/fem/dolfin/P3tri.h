// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.2.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __P3TRI_H
#define __P3TRI_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>

namespace dolfin
{
  
  class P3tri : public dolfin::FiniteElement
  {
  public:
  
    P3tri() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~P3tri()
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
      static unsigned int edge_reordering_0[2][2] = {{0, 1}, {1, 0}};
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
      int alignment = cell.edgeAlignment(0);
      int offset = mesh.numVertices();
      nodes[3] = offset + 2*cell.edgeID(0) + edge_reordering_0[alignment][0];
      nodes[4] = offset + 2*cell.edgeID(0) + edge_reordering_0[alignment][1];
      alignment = cell.edgeAlignment(1);
      nodes[5] = offset + 2*cell.edgeID(1) + edge_reordering_0[alignment][0];
      nodes[6] = offset + 2*cell.edgeID(1) + edge_reordering_0[alignment][1];
      alignment = cell.edgeAlignment(2);
      nodes[7] = offset + 2*cell.edgeID(2) + edge_reordering_0[alignment][0];
      nodes[8] = offset + 2*cell.edgeID(2) + edge_reordering_0[alignment][1];
      offset = offset + 2*mesh.numEdges();
      nodes[9] = offset + cell.id();
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(6.666666666666667e-01, 3.333333333333333e-01);
      points[4] = map(3.333333333333334e-01, 6.666666666666666e-01);
      points[5] = map(0.000000000000000e+00, 6.666666666666667e-01);
      points[6] = map(0.000000000000000e+00, 3.333333333333334e-01);
      points[7] = map(3.333333333333333e-01, 0.000000000000000e+00);
      points[8] = map(6.666666666666666e-01, 0.000000000000000e+00);
      points[9] = map(3.333333333333333e-01, 3.333333333333333e-01);
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
      FiniteElementSpec s("Lagrange", "triangle", 3);
      return s;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
  };
  
}

#endif
