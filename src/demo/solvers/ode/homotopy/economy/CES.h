// Copyright (C) 2005 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#ifndef __CES_H
#define __CES_H

#include <dolfin.h>
#include "Economy.h"

using std::pow;

using namespace dolfin;

// Constant elasticity of substitution (CES) economy (rational form)

class CES : public Economy
{
public:

  CES(unsigned int m, unsigned int n, real epsilon) : Economy(m, n)
  {
    // Choose b such that epsilon < b_i < 1 for all i
    b = new real[m];
    for (unsigned int i = 0; i < m; i++)
      b[i] = epsilon + dolfin::rand()*(1.0 - epsilon);

    /*
      a[0][0] = 2.0; a[0][1] = 1.0;
      a[1][0] = 0.5; a[1][1] = 1.0;
      
      w[0][0] = 1.0; w[0][1] = 0.0;
      w[1][0] = 0.0; w[1][1] = 1.0;
      
      b[0] = 0.1;
      b[1] = 0.1;
    */

    init(&tmp0);
    init(&tmp1);
  }

  ~CES() { delete [] b; }

  /*
  complex z0(unsigned int i)
  {
    if ( i == 0 )
      return 1.01362052581566e-02;
    else
      return 9.89863794741843e-01;
  }
  */
  
  void F(const complex z[], complex y[])
  {
    // First equation: normalization
    y[0] = sum(z) - 1.0;
    
    // Precompute scalar products
    for (unsigned int i = 0; i < m; i++)
      tmp0[i] = dot(w[i], z) / bdot(a[i], z, 1.0 - b[i]);
    
    // Evaluate right-hand side
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
	sum += a[i][j] * pow(z[j], -b[i]) * tmp0[i] - w[i][j];
      y[j] = sum;
    }
  }

  void JF(const complex z[], const complex x[], complex y[])
  {
    // First equation: normalization
    y[0] = sum(x);

    // First term
    for (unsigned int i = 0; i < m; i++)
    {
      const complex wx = dot(w[i], x);
      const complex az = bdot(a[i], z, 1.0 - b[i]);
      tmp0[i] = wx / az;
    }
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
	sum += a[i][j] * pow(z[j], -b[i]) * tmp0[i];
      y[j] = sum;
    }

    // Second term
    for (unsigned int i = 0; i < m; i++)
    {
      const complex wz = dot(w[i], z);
      const complex az = bdot(a[i], z, 1.0 - b[i]);
      tmp0[i] = b[i] * wz / az;
    }
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
	sum += a[i][j] * pow(z[j], -1.0 - b[i]) * x[j] * tmp0[i];
      y[j] -= sum;
    }

    // Third term
    for (unsigned int i = 0; i < m; i++)
    {
      const complex wz  = dot(w[i], z);
      const complex az  = bdot(a[i], z, 1.0 - b[i]);
      const complex axz = bdot(a[i], x, z, -b[i]);
      tmp0[i] = (1.0 - b[i]) * wz * axz / (az * az);
    }
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
	sum += a[i][j] * pow(z[j], -b[i]) * tmp0[i];
      y[j] -= sum;
    }
  }

  /*

  void G(const complex z[], complex y[])
  {
    // First equation: normalization
    y[0] = sum(z) - 1.0;
    
    // Only one consumer
    const unsigned int i = 0;

    // Precompute scalar product
    tmp0[i] = dot(w[i], z) / bdot(a[i], z, 1.0 - b[i]);
    
    // Evaluate right-hand side
    for (unsigned int j = 1; j < n; j++)
      y[j] = a[i][j] * pow(z[j], -b[i]) * tmp0[i] - w[i][j];
  }

  void JG(const complex z[], const complex x[], complex y[])
  {
    // First equation: normalization
    y[0] = sum(x);

    // Only one consumer
    const unsigned int i = 0;

    // First term
    complex wx = dot(w[i], x);
    complex az = bdot(a[i], z, 1.0 - b[i]);
    tmp0[i] = wx / az;
    for (unsigned int j = 1; j < n; j++)
      y[j] = a[i][j] * pow(z[j], -b[i]) * tmp0[i];

    // Second term
    complex wz = dot(w[i], z);
    az = bdot(a[i], z, 1.0 - b[i]);
    tmp0[i] = b[i] * wz / az;
    for (unsigned int j = 1; j < n; j++)
	y[j] -= a[i][j] * pow(z[j], -1.0 - b[i]) * x[j] * tmp0[i];
    
    // Third term
    wz  = dot(w[i], z);
    az  = bdot(a[i], z, 1.0 - b[i]);
    complex axz = bdot(a[i], x, z, -b[i]);
    tmp0[i] = (1.0 - b[i]) * wz * axz / (az * az);
    for (unsigned int j = 1; j < n; j++)
      y[j] = a[i][j] * pow(z[j], -b[i]) * tmp0[i];
  }

  */

  unsigned int degree(unsigned int i) const
  {
    if ( i == 0 )
      return 1;
    else
      return m;
  }
  
  // Vector of exponents
  real* b;
  
};

// Constant elasticity of substitution (CES) economy (polynomial form)

class PolynomialCES : public Economy
{
public:

  PolynomialCES(unsigned int m, unsigned int n, real epsilon) : Economy(m, n)
  {
    // Choose b such that epsilon < b_i < 1 for all i
    b = new real[m];
    for (unsigned int i = 0; i < m; i++)
      b[i] = epsilon + dolfin::rand()*(1.0 - epsilon);

    init(&tmp0);
    init(&tmp1);
    init(&tmp2);
    init(&tmp3);
  }

  ~PolynomialCES() { delete [] b; }

  void F(const complex z[], complex y[])
  {
    // First equation: normalization
    const complex zsum = sum(z);
    y[0] = zsum - 1.0;
    
    // Precompute scalar products
    real bsum = 0.0;
    for (unsigned int i = 0; i < m; i++)
    {
      tmp0[i] = dot(w[i], z);
      tmp1[i] = bdot(a[i], z, 1.0 - b[i]);
      bsum += b[i];
    }
    
    // Precompute product of all factors
    complex product = 1.0;
    for (unsigned int i = 0; i < m; i++)
      product *= tmp1[i];

    // Precompute dominating term
    complex extra = zsum;
    for (unsigned int i = 1; i < m; i++)
      extra *= zsum;
    extra -= 1.0;

    // Evaluate right-hand side
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      const complex tmp = pow(z[j], bsum) * product;
      for (unsigned int i = 0; i < m; i++)
      {
	const real di = bsum - b[i];
	sum += a[i][j] * tmp0[i] * pow(z[j], di) * product / tmp1[i];
	sum -= w[i][j] * tmp;
      }
      y[j] = sum + extra;
    }
  }

  void JF(const complex z[], const complex x[], complex y[])
  {
    // First equation: normalization
    const complex zsum = sum(z);
    const complex xsum = sum(x);
    y[0] = xsum;

    // Precompute scalar products
    real bsum = 0.0;
    for (unsigned int i = 0; i < m; i++)
    {
      tmp0[i] = dot(w[i], z);
      tmp1[i] = bdot(a[i], z, 1.0 - b[i]);
      tmp2[i] = dot(w[i], x);
      tmp3[i] = bdot(a[i], x, z, -b[i]);
      bsum += b[i];
    }
    
    // Precompute product of all factors
    complex product = 1.0;
    for (unsigned int i = 0; i < m; i++)
      product *= tmp1[i];

    // Precompute sum of all terms
    complex rsum = 0.0;
    for (unsigned int r = 0; r < m; r++)
      rsum += (1.0 - b[r]) * tmp3[r] * product / tmp1[r];

    // Precompute dominating term
    complex extra = static_cast<complex>(m) * xsum;
    for (unsigned int i = 1; i < m; i++)
      extra *= zsum;

    // Add terms of Jacobian
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
      {
	const real di = bsum - b[i];
	
	// First term
	sum += a[i][j]*(tmp2[i]*pow(z[j], di) + tmp0[i]*di*pow(z[j], di-1.0)*x[j]) *
	  product / tmp1[i];

	// Second term
	complex tmp = 0.0;
	for (unsigned int r = 0; r < m; r++)
	  if ( r != i )
	    tmp += (1.0 - b[r]) * tmp3[r] * product / (tmp1[r] * tmp1[i]);
	sum += a[i][j] * tmp0[i] * pow(z[j], di) * tmp;

	// Third term
	sum -= w[i][j] * bsum * pow(z[j], bsum - 1.0) * x[j] * product;

	// Fourth term
	sum -= w[i][j] * pow(z[j], bsum) * rsum;
      }
      y[j] = sum; // + extra;
    }
  }

  unsigned int degree(unsigned int i) const
  {
    if ( i == 0 )
      return 1;
    else
      return m;
  }
  
  // Vector of exponents
  real* b;

};

#endif
