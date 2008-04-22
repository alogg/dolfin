// Copyright (C) 2008 Solveig Bruvoll and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-03-28
// Last changed: 2008-04-22
//
// Transfinite mean value interpolation (ALE mesh smoothing)

#include <dolfin.h>
#include <string.h>
#include <cmath>

using namespace dolfin;

#define SIGN(v) (v < 0 ? -1 :1)
#define det(u,v,w) (u[0]*(v[1]*w[2]-v[2]*w[1])-u[1]*(v[0]*w[2]-v[2]*w[0])+u[2]*(v[0]*w[1]-v[1]*w[0]))
const real epsilon=1.0e-8;

unsigned int index(unsigned int i, unsigned int n) {
  if (i>=n) return i-n;
  if (i<0) return i+n;
  return i;
}

real dist(const real* x, const real* y, unsigned int dim)
{
  real s = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    s += (x[i] - y[i])*(x[i] - y[i]);
  return sqrt(s);
}

unsigned int next(unsigned int i, unsigned int dim)
{
  return (i == dim - 1 ? 0 : i + 1);
}

unsigned int previous(unsigned int i, unsigned int dim)
{
  return (i == 0 ? dim - 1 : i - 1);
}

void computeWeights(real* w, real** u, real* d, unsigned int dim, unsigned int num_vertices)
{
  real * ell = new real [num_vertices];
  real * theta = new real [num_vertices];
  real h=0;
  
  for (unsigned int i=0; i<num_vertices; i++) {
    const unsigned int ind1 = next(i, num_vertices);
    const unsigned int ind2 = previous(i, num_vertices);

    ell[i]=dist(u[ind1],u[ind2],dim); 
    		 
    theta[i]=2*asin(ell[i]/2.0);
    h+=theta[i]/2.0;
  }
  
    
  real sinus;
  real * c=new real[num_vertices];
  real * s = new real[num_vertices];
  
  for (unsigned int i=0; i<num_vertices; i++) {
    const unsigned int ind1 = next(i, num_vertices);
    const unsigned int ind2 = previous(i, num_vertices);

    c[i]=(2*sin(h)*sin(h-theta[i]))/(sin(theta[ind1])*sin(theta[ind2]))-1;
    sinus=1-c[i]*c[i];
    if(sinus<0 || sqrt(sinus)<epsilon) {
      for (unsigned int i=0; i<num_vertices; i++) 
	w[i]=0;
      return;
    }
    s[i]=SIGN(det(u[0],u[1],u[2]))*sqrt(sinus);
  }
  
  
  for (unsigned int i=0; i<num_vertices; i++) {
    const unsigned int ind1 = next(i, num_vertices);
    const unsigned int ind2 = previous(i, num_vertices);
    
    w[i]=(theta[i]-c[ind1]*theta[ind2]-c[ind2]*theta[ind1])/(d[i]*sin(theta[ind1])*s[ind2]);
  } 
  
}

void meanValue(real* new_x, unsigned int dim, Mesh& new_boundary,
	       Mesh& mesh, MeshFunction<unsigned int>& vertex_map,
	       Vertex& vertex)
{
  // Check if the point is on the boundary (no need to compute new coordinate)
  for (VertexIterator v(new_boundary); !v.end(); ++v)
  {
    if (vertex_map(*v) == vertex.index())
    {
      memcpy(new_x, v->x(), dim*sizeof(real));
      return;
    }
  }

  const real* old_x = vertex.x();
  cout << "Old x: " << old_x[0] << " " << old_x[1] <<" "<<old_x[2];
  
  const unsigned int size = new_boundary.numVertices();
  real * d = new real[size];
  real ** u = new real * [size];

  // Compute distance d and direction vector u from x to all p
  for (VertexIterator v(new_boundary);  !v.end(); ++v) {

    // Old position of point x
    const real* x = vertex.x();
    
    //Old position of vertex v in boundary:
    const real* p = mesh.geometry().x(vertex_map(*v));

    //distance from x to each point at the boundary
    d[v->index()] = dist(p, x, dim);
          
    //compute direction vector for p-x.    
    u[v->index()] = new real [dim];
    for (unsigned int i=0; i<dim; i++)
      u[v->index()][i]=(p[i] - x[i]) / d[v->index()];
        
  }
  const unsigned int num_vertices = new_boundary.topology().dim() + 1;
  
  // Local arrays
  real * w      = new real [num_vertices];
  real ** new_p = new real * [num_vertices];
  real * dCell  = new real [num_vertices];  
  real ** uCell = new real * [num_vertices];
  
  // Set new x to zero
  for (unsigned int i = 0; i < dim; ++i)
    new_x[i] = 0.0;
  
  // Iterate over all cells in boundary
  real totalW = 0.0;
  for (CellIterator c(new_boundary); !c.end(); ++c)
  {
    // Get local data
    unsigned int ind;
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      
      ind=v.pos();
      new_p[ind] = v->x();
      uCell[ind] = u[v->index()];
      dCell[ind] = d[v->index()];
    }
    
    //Compute weights w.
    computeWeights(w, uCell, dCell, dim, num_vertices);

    // Compute sum of weights
    for (unsigned int i=0; i<num_vertices; i++)
      totalW += w[i];
    
    // Compute new position
    for (unsigned int j=0; j<dim; j++)
      for (unsigned int i=0; i<num_vertices; i++)
	new_x[j] += w[i]*new_p[i][j];
      
    //psi[xnr]=1/totalW;
  }

  // Scale by totalW
  for (unsigned int i = 0; i < dim; i++)
    new_x[i] /= totalW;

  cout<<"  New x: "<<new_x[0]<<" "<<new_x[1]<<" "<<new_x[2]<<endl;

  // Free memory for d
  delete [] d;
  
  // Free memory for u
  for (unsigned int i = 0; i < size; ++i)
    delete [] u[i];
  delete [] u;
  
  // Free memory for local arrays
  delete [] w;
  delete [] new_p;
  delete [] uCell;
  delete [] dCell;
}

void deform(Mesh& mesh, Mesh& new_boundary, MeshFunction<unsigned int>& vertex_map, MeshFunction<unsigned int>& cell_map)
{
  // Extract old coordinates
  const unsigned int dim = mesh.geometry().dim();
  const unsigned int size = mesh.numVertices()*dim;
  real* new_x = new real[size];

  // Iterate over coordinates in mesh
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    meanValue(new_x + v->index()*dim, dim, new_boundary,
	      mesh, vertex_map, *v);
  }

  // Update mesh coordinates
  for (VertexIterator v(mesh); !v.end(); ++v)
    memcpy(v->x(), new_x + v->index()*dim, dim*sizeof(real));
  
  cout<<"mesh size: "<<mesh.numVertices()<<endl;
  cout<<"boundary size: "<<new_boundary.numVertices()<<endl;
  delete [] new_x;
}

int main()
{
  UnitCube mesh(3, 3, 3);
  plot(mesh);

  // For testing
  File mesh_file("mesh.xml");
  mesh_file << mesh;

  // Save mesh to file for visualization
  File f0("mesh_original.pvd");
  f0 << mesh;

  MeshFunction<unsigned int> vertex_map;
  MeshFunction<unsigned int> cell_map;
  BoundaryMesh boundary(mesh, vertex_map, cell_map);

  // For testing
  File boundary_file("boundary.xml");
  boundary_file << boundary;

  // Deform boundary
  for (VertexIterator v(boundary); !v.end(); ++v)
  {
    real* x = v->x();
    x[0] = x[0] + x[2];
  }
  plot(boundary);

  // Deform mesh
  deform(mesh, boundary, vertex_map, cell_map);
  plot(mesh);

  // Save mesh to file for visualization
  File f1("mesh_deformed.pvd");
  f1 << mesh;

  return 0;
}
