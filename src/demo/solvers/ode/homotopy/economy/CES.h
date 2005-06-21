// Copyright (C) 2005 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#ifndef __CES_H
#define __CES_H

#include <dolfin.h>
#include "Economy.h"

using std::pow;

using namespace dolfin;

/// Constant Elasticity of Substitution (CES) economy with m traders
/// and n goods. The system is implemented in three different
/// versions: rational form with rational exponents, polynomial form
/// with rational exponents, and polynomial form with integer
/// exponents (true polynomial form).

/// This class implements the original rational form with rational exponents.

class RationalRationalCES : public Economy
{
public:

  RationalRationalCES(unsigned int m, unsigned int n) : Economy(m, n)
  {
    b = new real[m];
    for (unsigned int i = 0; i < m; i++)
      b[i] = 0.0;
    
    init(&tmp0);
    init(&tmp1);
  }

  ~RationalRationalCES() { delete [] b; }

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

/// This class implements the polynomial form with rational exponents.

class PolynomialRationalCES : public Economy
{
public:

  PolynomialRationalCES(unsigned int m, unsigned int n) : Economy(m, n)
  {
    b = new real[m];
    for (unsigned int i = 0; i < m; i++)
      b[i] = 0.0;

    init(&tmp0);
    init(&tmp1);
    init(&tmp2);
    init(&tmp3);
  }

  ~PolynomialRationalCES() { delete [] b; }

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
    //complex extra = zsum;
    //for (unsigned int i = 1; i < m; i++)
    //  extra *= zsum;
    //extra -= 1.0;

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
      y[j] = sum; // + extra;
    }
  }

  void JF(const complex z[], const complex x[], complex y[])
  {
    // First equation: normalization
    //const complex zsum = sum(z);
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
    //complex extra = static_cast<complex>(m) * xsum;
    //for (unsigned int i = 1; i < m; i++)
    //  extra *= zsum;

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

/// This class implements the polynomial form with integer exponents.
/// Note that the vector beta is used to represent the substituted
/// values beta_i = b_i / epsilon and alpha = 1 / epsilon, with beta_i
/// and alpha integers.

class PolynomialIntegerCES : public Economy
{
public:

  PolynomialIntegerCES(unsigned int m, unsigned int n, bool real_valued = false) : Economy(m, n)
  {
    alpha = 1;
    beta = new unsigned int[m];
    for (unsigned int i = 0; i < m; i++)
    {
      beta[i] = 1;
    }
    
    this->real_valued = real_valued;

    init(&tmp0);
    init(&tmp1);
    init(&tmp2);
    init(&tmp3);
  }

  ~PolynomialIntegerCES() { delete [] beta; }

  // System F(z) = 0
  void F(const complex z[], complex y[])
  {
    // First equation: normalization
    const complex zsum = bsum(z, alpha);
    y[0] = zsum - 1.0;
    
    // Precompute scalar products
    unsigned int bsum = 0;
    for (unsigned int i = 0; i < m; i++)
    {
      tmp0[i] = bdot(w[i], z, alpha);
      tmp1[i] = bdot(a[i], z, alpha - beta[i]);
      bsum += beta[i];
    }
    
    // Precompute product of all factors
    complex product = 1.0;
    for (unsigned int i = 0; i < m; i++)
      product *= tmp1[i];

    // Evaluate right-hand side
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      const complex tmp = pow(z[j], bsum) * product;
      for (unsigned int i = 0; i < m; i++)
      {
	sum += a[i][j] * tmp0[i] * pow(z[j], bsum - beta[i]) * product / tmp1[i];
	sum -= w[i][j] * tmp;
      }
      y[j] = sum;
    }
  }

  // Jacobian dF/dz of system F(z) = 0
  void JF(const complex z[], const complex x[], complex y[])
  {
    // First equation: normalization
    const complex xsum = static_cast<real>(alpha) * bdot(x, z, alpha - 1);
    y[0] = xsum;

    // Precompute scalar products
    unsigned int bsum = 0;
    for (unsigned int i = 0; i < m; i++)
    {
      tmp0[i] = bdot(w[i], z, alpha);
      tmp1[i] = bdot(a[i], z, alpha - beta[i]);
      tmp2[i] = bdot(w[i], x, z, alpha - 1);
      tmp3[i] = bdot(a[i], x, z, alpha - beta[i] - 1);
      bsum += beta[i];
    }
    
    // Precompute product of all factors
    complex product = 1.0;
    for (unsigned int i = 0; i < m; i++)
      product *= tmp1[i];

    // Precompute sum of all terms
    complex rsum = 0.0;
    for (unsigned int r = 0; r < m; r++)
      rsum += static_cast<real>(alpha - beta[r]) * tmp3[r] * product / tmp1[r];

    // Add terms of Jacobian
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
      {
	// First term and second terms
	sum += a[i][j]*(static_cast<real>(alpha)*tmp2[i]*pow(z[j], bsum - beta[i]) + 
			tmp0[i]*static_cast<real>(bsum - beta[i])*
			pow(z[j], bsum - beta[i] - 1)*x[j]) * product / tmp1[i];

	// Third term
	complex tmp = 0.0;
	for (unsigned int r = 0; r < m; r++)
	  if ( r != i )
	    tmp += static_cast<real>(alpha - beta[r]) * tmp3[r] * product / (tmp1[r] * tmp1[i]);
	sum += a[i][j] * tmp0[i] * pow(z[j], bsum - beta[i]) * tmp;

	// Forth term
	sum -= w[i][j] * static_cast<real>(bsum) * pow(z[j], bsum - 1) * x[j] * product;

	// Fifth term
	sum -= w[i][j] * pow(z[j], bsum) * rsum;
      }
      y[j] = sum;
    }
  }

  void modify(complex z[])
  {
    for (unsigned int j = 0; j < n; j++)
    {
      // Scale back
      z[j] = std::pow(z[j], alpha);
      
      // Set almost zero imaginary parts to zero
      if ( std::abs(z[j].imag()) < DOLFIN_EPS )
	z[j] = z[j].real();
    }
  }

  bool verify(const complex z[])
  {
    const real tol = 2e-12;
    bool ok = true;

    dolfin_info("Verifying solution:");

    // Check normalization
    if ( std::abs(sum(z) - 1.0) < tol )
      dolfin_info("  - Normalization:  ok");
    else
    {
      ok = false;
      dolfin_info("  - Normalization:  failed");
    }

    // Precompute scalar products
    for (unsigned int i = 0; i < m; i++)
    {
      const real bi = static_cast<real>(beta[i]) / static_cast<real>(alpha);
      tmp0[i] = dot(w[i], z) / bdot(a[i], z, 1.0 - bi);
    }    

    // Check first equation (the one replaced with normalization)
    complex sum = 0.0;
    for (unsigned int i = 0; i < m; i++)
    {
      const real bi = static_cast<real>(beta[i]) / static_cast<real>(alpha);
      sum += a[i][0] * tmp0[i] / pow(z[0], bi) - w[i][0];
    }
    if ( std::abs(sum) < tol )
      dolfin_info("  - First equation: ok");
    else
    {
      ok = false;
      dolfin_info("  - First equation: failed");
    }

    // Check remaining equations
    real maxsum = 0.0;
    for (unsigned int j = 1; j < n; j++)
    {
      complex sum = 0.0;
      for (unsigned int i = 0; i < m; i++)
      {
	const real bi = static_cast<real>(beta[i]) / static_cast<real>(alpha);
	sum += a[i][j] * tmp0[i] / pow(z[j], bi) - w[i][j];
      }
      maxsum = std::max(maxsum, std::abs(sum));
    }
    if ( std::abs(sum) < tol )
      dolfin_info("  - Rest of system: ok");
    else
    {
      ok = false;
      dolfin_info("  - Rest of system: failed");
    }

    // Check if solution is real-valued
    if ( real_valued )
    {
      bool all_real = true;
      for (unsigned int j = 0; j < n; j++)
      {
	if ( std::abs(z[j].imag()) > tol )
	{
	  all_real = false;
	  break;
	}
      }
      if ( all_real )
	dolfin_info("  - Real-valued:    ok");
      else
      {
	ok = false;
	dolfin_info("  - Real-valued:    failed");
      }
    }

    return ok;
  }

  unsigned int degree(unsigned int i) const
  {
    if ( i == 0 )
      return 1;
    else
      return m;
  }
  
  // Scaled exponents (substituted integer values)
  unsigned int* beta;

  // Scaled exponent of numerator (1 / epsilon)
  unsigned int alpha;

private:

  // True if we only want real-valued solutions
  bool real_valued;

};

#endif
