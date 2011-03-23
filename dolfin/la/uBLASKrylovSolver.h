// Copyright (C) 2006-2009 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2006-2010.
//
// First added:  2006-05-31
// Last changed: 2011-03-24

#ifndef __UBLAS_KRYLOV_SOLVER_H
#define __UBLAS_KRYLOV_SOLVER_H

#include <set>
#include <string>
#include <boost/shared_ptr.hpp>
#include <dolfin/common/types.h>
#include "ublas.h"
#include "GenericLinearSolver.h"
#include "uBLASKrylovMatrix.h"
#include "uBLASMatrix.h"
#include "uBLASVector.h"
#include "uBLASPreconditioner.h"

namespace dolfin
{

  class GenericMatrix;
  class GenericVector;

  /// This class implements Krylov methods for linear systems
  /// of the form Ax = b using uBLAS data types.

  class uBLASKrylovSolver : public GenericLinearSolver
  {
  public:

    /// Create Krylov solver for a particular method and preconditioner
    uBLASKrylovSolver(std::string solver_type="default",
                      std::string pc_type="default");

    /// Create Krylov solver for a particular uBLASPreconditioner
    uBLASKrylovSolver(uBLASPreconditioner& pc);

    /// Create Krylov solver for a particular method and uBLASPreconditioner
    uBLASKrylovSolver(std::string solver_type,
                      uBLASPreconditioner& preconditioner);

    /// Destructor
    ~uBLASKrylovSolver();

    /// Solve the operator (matrix)
    void set_operator(const GenericMatrix& A)
    { error("set_operator(A) is not implemented."); }

    /// Return the operator (matrix)
    const GenericMatrix& get_operator() const
    {
      error("get_operator() is not implemented.");
      return *(static_cast<GenericMatrix*>(0)); // code will not be reached
    }

    /// Solve linear system Ax = b and return number of iterations
    uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b);

    /// Solve linear system Ax = b and return number of iterations (dense matrix)
    uint solve(const uBLASMatrix<ublas_dense_matrix>& A, uBLASVector& x,
               const uBLASVector& b);

    /// Solve linear system Ax = b and return number of iterations (sparse matrix)
    uint solve(const uBLASMatrix<ublas_sparse_matrix>& A, uBLASVector& x,
               const uBLASVector& b);

    /// Solve linear system Ax = b and return number of iterations (virtual matrix)
    uint solve(const uBLASKrylovMatrix& A, uBLASVector& x, const uBLASVector& b);

    /// Default parameter values
    static Parameters default_parameters();

  private:

    /// Select solver and solve linear system Ax = b and return number of iterations
    template<class Mat>
    uint solve_krylov(const Mat& A, uBLASVector& x, const uBLASVector& b);

    /// Solve linear system Ax = b using CG
    template<class Mat>
    uint solveCG(const Mat& A, uBLASVector& x, const uBLASVector& b,
                 bool& converged) const;

    /// Solve linear system Ax = b using restarted GMRES
    template<class Mat>
    uint solveGMRES(const Mat& A, uBLASVector& x, const uBLASVector& b,
                        bool& converged) const;

    /// Solve linear system Ax = b using BiCGStab
    template<class Mat>
    uint solveBiCGStab(const Mat& A, uBLASVector& x, const uBLASVector& b,
                        bool& converged) const;

    /// Select and create named preconditioner
    void select_preconditioner(std::string pc_type);

    /// Read solver parameters
    void read_parameters();

    /// Krylov method
    std::string solver_type;

    /// Preconditioner
    boost::shared_ptr<uBLASPreconditioner> pc;

    // Available solver types
    static const std::set<std::string> solver_types;

    /// Solver parameters
    double rtol, atol, div_tol;
    uint max_it, restart;
    bool report;

    /// True if we have read parameters
    bool parameters_read;

  };
  //---------------------------------------------------------------------------
  // Implementation of template functions
  //---------------------------------------------------------------------------
  template<class Mat>
  dolfin::uint uBLASKrylovSolver::solve_krylov(const Mat& A,
                                               uBLASVector& x,
                                               const uBLASVector& b)
  {
    if (solver_types.count(solver_type) == 0)
      error("Requested Krylov solver '%s' not available in uBLASKrylovSolver.", solver_type.c_str());

    // Check dimensions
    uint M = A.size(0);
    uint N = A.size(1);
    if ( N != b.size() )
      error("Non-matching dimensions for linear system.");

    // Reinitialise x if necessary
    if (x.size() != b.size())
    {
      x.resize(b.size());
      x.zero();
    }

    // Read parameters if not done
    if ( !parameters_read )
      read_parameters();

    // Write a message
    if ( report )
      info("Solving linear system of size %d x %d (uBLAS Krylov solver).", M, N);

    // Initialise preconditioner if necessary
    pc->init(A);

    // Choose solver and solve
    bool converged = false;
    uint iterations = 0;
    if (solver_type == "cg")
      iterations = solveCG(A, x, b, converged);
    else if (solver_type == "gmres")
      iterations = solveGMRES(A, x, b, converged);
    else if (solver_type == "bicgstab")
      iterations = solveBiCGStab(A, x, b, converged);
    else if (solver_type == "default")
      iterations = solveBiCGStab(A, x, b, converged);
    else
      error("Requested Krylov method unknown.");

    // Check for convergence
    if (!converged)
    {
      bool error_on_nonconvergence = parameters["error_on_nonconvergence"];
      if (error_on_nonconvergence)
        error("uBLAS Krylov solver failed to converge.");
      else
        warning("uBLAS Krylov solver failed to converge.");
    }
    else if (report)
      info("Krylov solver converged in %d iterations.", iterations);

    return iterations;
  }
  //-----------------------------------------------------------------------------
  template<class Mat>
  dolfin::uint uBLASKrylovSolver::solveCG(const Mat& A,
                                          uBLASVector& x,
                                          const uBLASVector& b,
                                          bool& converged) const
  {
    warning("Conjugate-gradient method not yet programmed for uBLASKrylovSolver. Using GMRES.");
    return solveGMRES(A, x, b, converged);
  }
  //-----------------------------------------------------------------------------
  template<class Mat>
  dolfin::uint uBLASKrylovSolver::solveGMRES(const Mat& A, uBLASVector& x,
                                             const uBLASVector& b, bool& converged) const
  {
    // Get underlying uBLAS vectors
    ublas_vector& _x = x.vec();
    const ublas_vector& _b = b.vec();

    // Get size of system
    const uint size = A.size(0);

    // Create residual vector
    uBLASVector r(size);
    ublas_vector& _r = r.vec();

    // Create H matrix and h vector
    ublas_matrix_cmajor_tri H(restart, restart);
    ublas_vector _h(restart+1);

    // Create gamma vector
    ublas_vector _gamma(restart+1);

    // Matrix containing v_k as columns.
    ublas_matrix_cmajor V(size, restart+1);

    // w vector
    uBLASVector w(size);
    ublas_vector& _w = w.vec();

    // Givens vectors
    ublas_vector _c(restart), _s(restart);

    // Miscellaneous storage
    double nu, temp1, temp2, r_norm = 0.0, beta0 = 0;

    converged = false;
    uint iteration = 0;
    while (iteration < max_it && !converged)
    {
      // Compute residual r = b -A*x
      //noalias(r) = b;
      //axpy_prod(A, -x, r, false);
      A.mult(x, r);
      _r *= -1.0;
      noalias(_r) += _b;

      // Apply preconditioner (use w for temporary storage)
      _w.assign(_r);
      pc->solve(r, w);

      // L2 norm of residual (for most recent restart)
      const double beta = norm_2(_r);

     // Save intial residual (from restart 0)
     if(iteration == 0)
       beta0 = beta;

     if(beta < atol)
     {
       converged = true;
       return iteration;
     }

      // Intialise gamma
      _gamma.clear();
      _gamma(0) = beta;

      // Create first column of V
      noalias(column(V, 0)) = _r/beta;

      // Modified Gram-Schmidt procedure
      uint subiteration = 0;
      uint j = 0;
      while (subiteration < restart && iteration < max_it && !converged && r_norm/beta < div_tol)
      {
        // Compute product w = A*V_j (use r for temporary storage)
        //axpy_prod(A, column(V, j), w, true);
        noalias(_r) = column(V, j);
        A.mult(r, w);

        // Apply preconditioner (use r for temporary storage)
        _r.assign(_w);
        pc->solve(w, r);

        for (uint i=0; i <= j; ++i)
        {
          _h(i)= inner_prod(_w, column(V,i));
          noalias(_w) -= _h(i)*column(V,i);
        }
        _h(j+1) = norm_2(_w);

        // Insert column of V (inserting v_(j+1)
        noalias(column(V,j+1)) = _w/_h(j+1);

        // Apply previous Givens rotations to the "new" column
        // (this could be improved? - use more uBLAS functions.
        //  The below has been taken from old DOLFIN code.)
        for(uint i=0; i<j; ++i)
        {
          temp1 = _h(i);
          temp2 = _h(i+1);
          _h(i)   = _c(i)*temp1 - _s(i)*temp2;
          _h(i+1) = _s(i)*temp1 + _c(i)*temp2 ;
        }

        // Compute new c_i and s_i
        nu = sqrt( _h(j)*_h(j) + _h(j+1)*_h(j+1) );

        // Direct access to c & s below leads to some strange compiler errors
        // when using vector expressions and noalias(). By using "subrange",
        // we are working with vector expressions rather than reals
        //c(j) =  h(j)/nu;
        //s(j) = -h(j+1)/nu;
        subrange(_c, j,j+1) =  subrange(_h, j,j+1)/nu;
        subrange(_s, j,j+1) = -subrange(_h, j+1,j+2)/nu;

        // Apply new rotation to last column
        _h(j)   = _c(j)*_h(j) - _s(j)*_h(j+1);
        _h(j+1) = 0.0;

        // Apply rotations to gamma
        temp1 = _c(j)*_gamma(j) - _s(j)*_gamma(j+1);
        _gamma(j+1) = _s(j)*_gamma(j) + _c(j)*_gamma(j+1);
        _gamma(j) = temp1;
        r_norm = fabs(_gamma(j+1));

        // Add h to H matrix. Would ne nice to use
        //   noalias(column(H, j)) = subrange(h, 0, restart);
        // but this gives an error when uBLAS debugging is turned onand H
        // is a triangular matrix
        for(uint i=0; i<j+1; ++i)
          H(i,j) = _h(i);

        // Check for convergence
        if( r_norm/beta0 < rtol || r_norm < atol )
          converged = true;

        ++iteration;
        ++subiteration;
        ++j;
      }

      // Eliminate extra rows and columns (this does not resize or copy, just addresses a range)
      ublas_matrix_cmajor_tri_range Htrunc(H, ublas::range(0,subiteration), ublas::range(0,subiteration) );
      ublas_vector_range _g(_gamma, ublas::range(0,subiteration));

      // Solve triangular system H*g and return result in g
      ublas::inplace_solve(Htrunc, _g, ublas::upper_tag ());

      // x_m = x_0 + V*y
      ublas_matrix_cmajor_range _v( V, ublas::range(0,V.size1()), ublas::range(0,subiteration) );
      axpy_prod(_v, _g, _x, false);
    }
    return iteration;
  }
  //-----------------------------------------------------------------------------
  template<class Mat>
  dolfin::uint uBLASKrylovSolver::solveBiCGStab(const Mat& A, uBLASVector& x,
                                                const uBLASVector& b,
                                                bool& converged) const
  {
    // Get uderlying uBLAS vectors
    ublas_vector& _x = x.vec();
    const ublas_vector& _b = b.vec();

   // Get size of system
    const uint size = A.size(0);

    // Allocate vectors
    uBLASVector r(size), rstar(size), p(size), s(size), v(size), t(size), y(size), z(size);
    ublas_vector& _r = r.vec();
    ublas_vector& _rstar = rstar.vec();
    ublas_vector& _p = p.vec();
    ublas_vector& _s = s.vec();
    ublas_vector& _v = v.vec();
    ublas_vector& _t = t.vec();
    ublas_vector& _y = y.vec();
    ublas_vector& _z = z.vec();

    double alpha = 1.0, beta = 0.0, omega = 1.0, r_norm = 0.0;
    double rho_old = 1.0, rho = 1.0;

    // Compute residual r = b -A*x
    //r.assign(b);
    //axpy_prod(A, -x, r, false);
    A.mult(x, r);
    r *= -1.0;
    noalias(_r) += _b;

    const double r0_norm = norm_2(_r);
    if( r0_norm < atol )
    {
      converged = true;
      return 0;
    }

    // Initialise r^star, v and p
    _rstar.assign(_r);
    _v.clear();
    _p.clear();

    // Apply preconditioner to r^start. This is a trick to avoid problems in which
    // (r^start, r) = 0  after the first iteration (such as PDE's with homogeneous
    // Neumann bc's and no forcing/source term.
    pc->solve(rstar, r);

    // Right-preconditioned Bi-CGSTAB

    // Start iterations
    converged = false;
    uint iteration = 0;
    while (iteration < max_it && !converged && r_norm/r0_norm < div_tol)
    {
      // Set rho_n = rho_n+1
      rho_old = rho;

      // Compute new rho
      rho = ublas::inner_prod(_r, _rstar);
      if( fabs(rho) < 1e-25 )
        error("BiCGStab breakdown. rho = %g", rho);

      beta = (rho/rho_old)*(alpha/omega);

      // p = r1 + beta*p - beta*omega*A*p
      p *= beta;
      noalias(_p) += _r - beta*omega*_v;

      // My = p
      pc->solve(y, p);

      // v = A*y
      //axpy_prod(A, y, v, true);
      A.mult(y, v);

      // alpha = (r, rstart) / (v, rstar)
      alpha = rho/ublas::inner_prod(_v, _rstar);

      // s = r - alpha*v
      noalias(_s) = _r - alpha*_v;

      // Mz = s
      pc->solve(z, s);

      // t = A*z
      //axpy_prod(A, z, t, true);
      A.mult(z, t);

      // omega = (t, s) / (t,t)
      omega = ublas::inner_prod(_t, _s)/ublas::inner_prod(_t, _t);

      // x = x + alpha*p + omega*s
      noalias(_x) += alpha*_y + omega*_z;

      // r = s - omega*t
      noalias(_r) = _s - omega*_t;

      // Compute norm of the residual and check for convergence
      r_norm = norm_2(_r);
      if( r_norm/r0_norm < rtol || r_norm < atol)
        converged = true;

      ++iteration;
    }

    return iteration;
  }
  //-----------------------------------------------------------------------------
}

#endif
