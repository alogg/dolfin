// Copyright (C) 2005 Johan Jansson.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2004-2005.
//
// First added:  2005
// Last changed: 2005

//#include <iostream>
#include <sstream>
#include <iomanip>

#include "dolfin/timeinfo.h"
#include "dolfin/ElasticityUpdatedSolver.h"
#include "dolfin/ElasticityUpdated.h"
#include "dolfin/ElasticityUpdatedSigma0.h"
#include "dolfin/ElasticityUpdatedSigma1.h"
#include "dolfin/ElasticityUpdatedSigma2.h"
#include "dolfin/ElasticityUpdatedProj.h"
#include "dolfin/ElasticityUpdatedMass.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
ElasticityUpdatedSolver::ElasticityUpdatedSolver(Mesh& mesh, 
						 Function& f,
						 Function& v0,
						 Function& rho,
						 real E, real nu, real nuv,
						 real nuplast,
						 BoundaryCondition& bc,
						 real k, real T)
  : mesh(mesh), f(f), v0(v0), rho(rho), E(E), nu(nu), nuv(nuv),
    nuplast(nuplast), bc(bc), k(k),
    T(T), counter(0), lastsample(0),
    lambda(E * nu / ((1 + nu) * (1 - 2 * nu))),
    mu(E / (2 * (1 + nu))),
    t(0.0), rtol(1.0e-4), maxiters(10), do_plasticity(false), yield(0.0),
    v1(x2_1, mesh, element1),
    sigma0_1(xsigma0_1, mesh, element2_0),
    sigma1_1(xsigma1_1, mesh, element2_1),
    sigma2_1(xsigma2_1, mesh, element2_2),
    sigma0_0(xsigma0_0, mesh, element2_0),
    sigma1_0(xsigma1_0, mesh, element2_1),
    sigma2_0(xsigma2_0, mesh, element2_2),
    epsilon0_1(xepsilon0_1, mesh, element2_0),
    epsilon1_1(xepsilon1_1, mesh, element2_1),
    epsilon2_1(xepsilon2_1, mesh, element2_2),
    sigmanorm(xsigmanorm, mesh, element2_2),
    Lv(f, sigma0_1, sigma1_1, sigma2_1,
       epsilon0_1, epsilon1_1, epsilon2_1, nuv),
    Lsigma0(v1, sigma0_1, sigma1_1, sigma2_1,
	    sigmanorm, lambda, mu, nuplast),
    Lsigma1(v1, sigma0_1, sigma1_1, sigma2_1,
	    sigmanorm, lambda, mu, nuplast),
    Lsigma2(v1, sigma0_1, sigma1_1, sigma2_1,
	    sigmanorm, lambda, mu, nuplast)

  
{
  // Do nothing
  init();
}
//-----------------------------------------------------------------------------
void ElasticityUpdatedSolver::init()
{
  Matrix M;

  // FIXME: Temporary fix, sizes should be automatically computable
  int Nv = FEM::size(mesh, element1);
  int Nsigma = FEM::size(mesh, element2_0);

  x1_0.init(Nv);
  x1_1.init(Nv);
  x2_0.init(Nv);
  x2_1.init(Nv);

  xtmp1.init(Nv);
  xtmp2.init(Nv);

  msigma.init(Nsigma);

  xtmp0_1.init(Nsigma);
  xtmp1_1.init(Nsigma);
  xtmp2_1.init(Nsigma);

  xsigma0_0.init(Nsigma);
  xsigma0_1.init(Nsigma);

  xsigma1_0.init(Nsigma);
  xsigma1_1.init(Nsigma);

  xsigma2_0.init(Nsigma);
  xsigma2_1.init(Nsigma);

  xepsilon0_1.init(Nsigma);
  xepsilon1_1.init(Nsigma);
  xepsilon2_1.init(Nsigma);

  xsigmanorm.init(Nsigma);

  xsigma0_1 = 0;
  xsigma1_1 = 0;
  xsigma2_1 = 0;

  stepresidual.init(Nv);

  // Set initial velocities

  AffineMap map;
  v0.set(element1);

  for(CellIterator c(&mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    // Use DOLFIN's interpolation

    real coefficients[element1.spacedim()];
    int dofs[element1.spacedim()];

    map.update(cell);
    v0.interpolate(coefficients, map);
    element1.dofmap(dofs, cell, mesh);

    for(uint i = 0; i < element1.spacedim(); i++)
      x2_1(dofs[i]) = coefficients[i];
  }

  cout << "x2_1:" << endl;
  x2_1.disp();

  xsigmanorm = 1.0;

  ElasticityUpdatedMass::BilinearForm amass(rho);

  // Assemble mass matrix
  FEM::assemble(amass, M, mesh);

  // Lump mass matrix
  FEM::lump(M, m);

  // Compute mass vector (sigma)
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    Cell& cell = *c;

    int dofs[element2_0.spacedim()];
    element2_0.dofmap(dofs, cell, mesh);

    real factor = 1.0 / cell.volume(); 

    for(uint i = 0; i < element2_0.spacedim(); i++)
      msigma(dofs[i]) = factor;
  }
}
//-----------------------------------------------------------------------------
void ElasticityUpdatedSolver::preparestep()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void ElasticityUpdatedSolver::step()
{
  // Make time step
  x1_0 = x1_1;
  x2_0 = x2_1;
  xsigma0_0 = xsigma0_1;
  xsigma1_0 = xsigma1_1;
  xsigma2_0 = xsigma2_1;
  
  t += k;
  f.set(t);
  //     cout << "t: " << t << endl;
  
  for(int iter = 0; iter < maxiters; iter++)
  {
    // Compute norm of stress (sigmanorm)
    if(do_plasticity)
    {
      for (CellIterator c(mesh); !c.end(); ++c)
      {
	Cell& cell = *c;

	int dofs[element2_0.spacedim()];
	element2_0.dofmap(dofs, cell, mesh);

	real proj = 1;
	real norm = 0;
	for(uint i = 0; i < element2_0.spacedim(); i++)
	{
	  norm = std::max(norm, fabs(xsigma0_0(dofs[i])));
	  norm = std::max(norm, fabs(xsigma1_0(dofs[i])));
	  norm = std::max(norm, fabs(xsigma2_0(dofs[i])));
	}
	
	if(norm > yield)
	{
	  cout << "sigmanorm(" << cell.id() << "): " << norm << endl;
	  proj = 1.0 / norm;
	}
	
	xsigmanorm(dofs[0]) = proj;
      }
    }
    
    //       dolfin_debug("Assembling sigma vectors");
    //       tic();
    
    // Assemble sigma0 vectors
    FEM::assemble(Lsigma0, xsigma0_1, mesh);
    FEM::assemble(Lsigma1, xsigma1_1, mesh);
    FEM::assemble(Lsigma2, xsigma2_1, mesh);
    
    VecPointwiseMult(xtmp0_1.vec(), xsigma0_1.vec(), msigma.vec());
    VecPointwiseMult(xtmp1_1.vec(), xsigma1_1.vec(), msigma.vec());
    VecPointwiseMult(xtmp2_1.vec(), xsigma2_1.vec(), msigma.vec());
    
    xtmp0_1.apply();
    xtmp1_1.apply();
    xtmp2_1.apply();
    
    xsigma0_1 = xsigma0_0;
    xsigma0_1.axpy(k, xtmp0_1);
    
    xsigma1_1 = xsigma1_0;
    xsigma1_1.axpy(k, xtmp1_1);
    
    xsigma2_1 = xsigma2_0;
    xsigma2_1.axpy(k, xtmp2_1);
    
    xepsilon0_1 *= 1.0 / lambda;
    xepsilon1_1 *= 1.0 / lambda;
    xepsilon2_1 *= 1.0 / lambda;
    
    // Assemble v vector
    FEM::assemble(Lv, xtmp1, mesh);
    
    b = xtmp1;
    b *= k;
    
    VecPointwiseDivide(stepresidual.vec(), xtmp1.vec(), m.vec());
    stepresidual *= k;
    stepresidual.axpy(-1, x2_1);
    stepresidual.axpy(1, x2_0);
    
    
    x2_1 += stepresidual;
    
    x1_1 = x1_0;
    x1_1.axpy(k, x2_0);
    
    cout << "stepresidual: " << stepresidual.norm(Vector::linf) << endl;
    
    if(stepresidual.norm(Vector::linf) <= rtol && iter >= 0)
    {
      cout << "converged" << endl;
      break;
    }
    else if(iter == maxiters - 1)
    {
      cout << "did not converge" << endl;
    }
  }
  
  
  
  // Update the mesh
  for (NodeIterator n(&mesh); !n.end(); ++n)
  {
    Node& node = *n;

    node.coord().x += k * v1(node, 0);
    node.coord().y += k * v1(node, 1);
    node.coord().z += k * v1(node, 2);
  }
}
//-----------------------------------------------------------------------------
void ElasticityUpdatedSolver::solve()
{
  cout << "lambda: " << lambda << endl;
  cout << "mu: " << mu << endl;
  
  File         file("elasticity.m");

  // Save the solution
  save(mesh, file, t);

  // Start a progress session
  Progress p("Time-stepping");
  
  // Start time-stepping
  while ( true && t < T ) {
  
    preparestep();
    step();
    
    // Save the solution
    save(mesh, file, t);

    // Benchmark
//     FEM::assemble(Lsigma0, xsigma0_1, mesh);
    
    // Update progress
    p = t / T;
  }
}
//-----------------------------------------------------------------------------
void ElasticityUpdatedSolver::save(Mesh& mesh, File& solutionfile, real t)
{
  real samplefreq = 1.0 / 33.0;

  while(lastsample + samplefreq < t || t == 0.0)
  {
    std::ostringstream fileid, filename;
    fileid.fill('0');
    fileid.width(6);
    
    fileid << counter;
    
    filename << "mesh" << fileid.str() << ".xml.gz";
    
    cout << "writing: " << filename.str() << " at t: " << t << endl;
    
    std::string foo = filename.str();
    const char *fname = foo.c_str();
    
    File meshfile(fname);
    
    meshfile << mesh;
    counter++;

    cout << "lastsample: " << lastsample << " t: " << t << endl;

    lastsample = std::min(t, lastsample + samplefreq);

    if(t == 0.0)
    {
      break;
    }
  }

}
//-----------------------------------------------------------------------------
void ElasticityUpdatedSolver::solve(Mesh& mesh,
				    Function& f,
				    Function& v0,
				    Function& rho,
				    real E, real nu, real nuv,
				    real nuplast,
				    BoundaryCondition& bc,
				    real k, real T)
{
  ElasticityUpdatedSolver solver(mesh, f, v0, rho, E, nu, nuv, nuplast,
				 bc, k, T);
  solver.solve();
}
//-----------------------------------------------------------------------------
