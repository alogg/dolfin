// Copyright (C) 2004 Johan Hoffman.
// Licensed under the GNU GPL Version 2.

#ifndef __NSE_MOMENTUM_H
#define __NSE_MOMENTUM_H

#include <dolfin/PDE.h>

namespace dolfin {
  
  class NSE_Momentum : public PDE {
  public:
    
    NSE_Momentum(Function::Vector& source,
		 Function::Vector& uprevious,
		 Function& viscosity,
		 Function::Vector& convection,
		 Function& pressure) : PDE(3, 3)
      {
	add(f,  source);
	add(up, uprevious);
	add(nu,  diffusion);
	add(b,  convection);
	add(p,  pressure);

	C1 = 2.0;
	C2 = 1.0;
      }
    
    real lhs(ShapeFunction::Vector& u, ShapeFunction::Vector& v)
    {

      unorm = sqrt(sqr(up(0))+sqr(up(1))+sqr(up(2)));
      if ( (h/nu) > 1.0 ) d1 = C1 * (0.5 / sqrt( 1.0/sqr(k) + sqr(unorm/h) ));
      else d1 = C1 * sqr(h);

      if ( (h/nu) > 1.0 ) d2 = C2 * h;
      else d2 = C2 * sqr(h);
      
      return
	( (u,v)*(1.0/k) + 0.5 * 
	  (nu*((grad(u(0)),grad(v(0))) + (grad(u(1)),grad(v(1))) + (grad(u(2)),grad(v(2)))) + 
	   (b,grad(u(0)))*v(0) + (b,grad(u(1)))*v(1) + (b,grad(u(2)))*v(2) + 
	   d1*((b,grad(u(0)))*(b,grad(v(0))) + (b,grad(u(1)))*(b,grad(v(1))) + (b,grad(u(2)))*(b,grad(v(2)))) + 
	   d2*(u(0).dx()+u(1).dy()+u(2).dz())*(v(0).dx()+v(1).dy()+v(2).dz())) ) * dK

    }
    
    real rhs(ShapeFunction::Vector& v)
    {

      unorm = sqrt(sqr(up(0))+sqr(up(1))+sqr(up(2)));
      if ( (h/nu) > 1.0 ) d1 = C1 * (0.5 / sqrt( 1.0/sqr(k) + sqr(unorm/h) ));
      else d1 = C1 * sqr(h);

      if ( (h/nu) > 1.0 ) d2 = C2 * h;
      else d2 = C2 * sqr(h);
      
      return
	( (up,v)*(1.0/k) - 0.5 * 
	  (nu*((grad(up(0)),grad(v(0))) +  (grad(up(1)),grad(v(1))) +  (grad(up(2)),grad(v(2)))) +  
	   (b,grad(up(0)))*v(0) + (b,grad(up(1)))*v(1) + (b,grad(up(2)))*v(2) + 
	   d1*((b,grad(up(0)))*(b,grad(v(0))) + (b,grad(up(1)))*(b,grad(v(1))) + (b,grad(up(2)))*(b,grad(v(2)))) +   
	   d2*(up(0).dx()+up(1).dy()+up(2).dz())*(v(0).dx()+v(1).dy()+v(2).dz())) -
	  d1*(p.dx()*(b,grad(v(0))) + p.dy()*(b,grad(v(1))) + p.dz()*(b,grad(v(2)))) + 
	  p*(v(0).dx()+v(1).dy()+v(2).dz()) ) * dK

    }
    
  private:    
    ElementFunction::Vector f;   // Source term
    ElementFunction::Vector up;  // Velocity value at left end-point
    ElementFunction nu;          // Viscosity
    ElementFunction::Vector b;   // Convection = linearized velocity
    ElementFunction p;           // linearized pressure

    real d1,d2,C1,C2,unorm;
  };
  
}

#endif
