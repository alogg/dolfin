// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-06-03
// Last changed: 2006-10-23

#include <dolfin/common/constants.h>
#include <dolfin/common/real.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/math/Legendre.h>
#include "RadauQuadrature.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
RadauQuadrature::RadauQuadrature(unsigned int n) : GaussianQuadrature(n)
{
  init();

  if ( !check(2*n-2) )
    error("Radau quadrature not ok, check failed.");

  //message("Radau quadrature computed for n = %d, check passed.", n);
}
//-----------------------------------------------------------------------------
void RadauQuadrature::disp() const
{
  cout << "Radau quadrature points and weights on [-1,1] for n = " 
       << n << ":" << endl;

  cout << " i    points                   weights" << endl;
  cout << "-----------------------------------------------------" << endl;

  for (unsigned int i = 0; i < n; i++)
    message("%2d   %.16e   %.16e", i, to_double(points[i]), to_double(weights[i]));
}
//-----------------------------------------------------------------------------
void RadauQuadrature::computePoints()
{
  // Compute the Radau quadrature points in [-1,1] as -1 and the zeros
  // of ( Pn-1(x) + Pn(x) ) / (1+x) where Pn is the n:th Legendre
  // polynomial. Computation is a little different than for Gauss and
  // Lobatto quadrature, since we don't know of any good initial
  // approximation for the Newton iterations.
  
  // Special case n = 1
  if ( n == 1 )
  {
    points[0] = -1.0;
    return;
  }

  Legendre p(n);
  real x, dx, step, sign;
  
  // Set size of stepping for seeking starting points
  step = 1.0 / ( double(n-1) * 15.0 );
  
  // Set the first nodal point which is -1
  points[0] = -1.0;
  
  // Start at -1 + step
  x = -1.0 + step;
  
  // Set the sign at -1 + epsilon
  sign = ( (p.eval(n-1,x) + p(x)) > 0 ? 1.0 : -1.0 );
  
  // Compute the rest of the nodes by Newton's method
  for (unsigned int i = 1; i < n; i++) {
    
    // Step to a sign change
    while ( (p.eval(n-1, x) + p(x))*sign > 0.0 )
      x += step;
    
    // Newton's method
    do
    {
      dx = - (p.eval(n-1, x) + p(x)) / (p.ddx(n-1, x) + p.ddx(x));
      x  = x + dx;
    } while ( abs(dx) > real_epsilon() );
    
    // Set the node value
    points[i] = x;
    
    // Fix step so that it's not too large
    if ( step > (points[i] - points[i-1])/10.0 )
      step = (points[i] - points[i-1]) / 10.0;
    
    // Step forward
    sign = - sign;
    x += step;
    
  }
  
}
//-----------------------------------------------------------------------------

