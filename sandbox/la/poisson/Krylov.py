# -*- coding: utf-8 -*-

"""
@license: Copyright (C) 2005
Author: Åsmund Ødegård, Ola Skavhaug, Gunnar Staff

Simula Research Laboratory AS

This file is part of PyCC.

PyCC  free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

PyCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PyFDM; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


ConjGrad: A generic implementation of the (preconditioned) conjugate gradient method.
"""

"""
Solve Ax = b with the Krylov methods.

A: Matrix or discrete operator. Must support A*x to operate on x.

x: Anything that A can operate on, as well as support inner product
and multiplication/addition/subtraction with scalars and vectors. Also
incremental operators like __imul__ and __iadd__ might be needed

b: Right-hand side, probably the same type as x.
"""



#from debug import debug

def debug(string, level):
    print string


import numpy as _n
from math import sqrt

#debug("Deprecation warning: You have imported Krylov.py. Use Solvers.py instead",level=0)

def inner(u,v):
    """Compute innerproduct of u and v.
       It is not computed here, we just check type of operands and
       call a suitable innerproduct method"""

    if isinstance(u,_n.ndarray):
        # assume both are numarrays:
        return _n.dot(u,v)
    else:
        # assume that left operand implements inner
        return u.inner(v) 

"""Conjugate Gradient Methods"""

def conjgrad(A, x, b, tolerance=1.0E-05, relativeconv=False, maxit=500):
    """
    conjgrad(A,x,b): Solve Ax = b with the conjugate gradient method.

    @param A: See module documentation
    @param x: See module documentation
    @param b: See module documentation

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.
    """

    r = b - A*x
    p = 1.0*r
    r0 = inner(r,r)

    if relativeconv:
        tolerance *= sqrt(inner(b,b))
    iter = 0
    while sqrt(r0) > tolerance and iter <= maxit:
        w = A*p
        a = r0/inner(p,w)
        x += a*p
        r -= a*w
        r1 = inner(r,r)
        p *= r1/r0
        p += r
        r0 = r1
        iter += 1
    #debug("CG finished, iter: %d, ||e||=%e" % (iter,sqrt(r0)))
    return x

def precondconjgrad(B, A, x, b, tolerance=1.0E-05, relativeconv=False, maxiter=500):
    r"""Preconditioned Conjugate Gradients method. Algorithm described here; 
    http://www.math-linux.com/spip.php?article55"""
    r = b - A*x
    z = B*r
    d = 1.0*z

    rz = inner(r,z)

    if relativeconv: tolerance *= sqrt(rz)
    iter = 0
    while sqrt(rz) > tolerance and iter <= maxiter:
        z = A*d
        alpha = rz / inner(d,z)
        x += alpha*d
        r -= alpha*z
        z = B*r
        rz_prev = rz
        rz = inner(r,z)
        beta =  rz / rz_prev
        d = z + beta*d
        iter += 1
    return x



def old_precondconjgrad(B, A, x, b, tolerance=1.0E-05, relativeconv=False, maxiter=500):
    """
    conjgrad(B,A,x,b): Solve Ax = b with the preconditioned conjugate gradient method.

    @param B: Preconditioner supporting the __mul__ operator for operating on
    x.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """

    r0 = b - A*x
    r = 1.0*r0
    z = 1.0*r
    s = B*r
    p = 1.0*s
    rho = rho0 = inner(s, r)
    if relativeconv:
        tolerance *= sqrt(inner(B*b,b))
    iter = 0
    while sqrt(rho0) > tolerance and iter <= maxiter:
        z = A*p
        t = B*z
        g = inner(p,z)
        a = rho0/g
        x += a*p
        r -= a*z
        s -= a*t
        rho = inner(s, r)
        b = rho/rho0
        p *= b; p += s
        rho0 = rho
        iter += 1
    debug("PCG finished, iter: %d, ||rho||=%e" % (iter,sqrt(rho0)))
    return x


"""Bi-conjugate gradients method"""

def BiCGStab(A, x, b, tolerance=1.0E-05, relativeconv=False, maxiter=1000, info=False):
    """
    BiCGStab(A,x,b): Solve Ax = b with the Biconjugate gradient stabilized
    method.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.
    """

    r = b - A*x
    p = r.copy()
    rs= r.copy()
    rr = inner(r,rs)

    if relativeconv:
        tolerance *= sqrt(inner(b,b))
    iter = 0
    debug("r0=" +  str(sqrt(inner(r,r))), 0)
    while sqrt(inner(r,r)) > tolerance  and  iter < maxiter:
        Ap    = A*p
        alpha = rr/inner(rs,Ap)
        print "alpha ", alpha
        s     = r-alpha*Ap
        As    = A*s
        w     = inner(As,s)/inner(As,As)
        x    += alpha*p+w*s
        r     = s -w*As
#        r = b - A*x
        rrn   = inner(r,rs)
        beta  = (rrn/rr)*(alpha/w)
        if beta==0.0:
            debug("BiCGStab breakdown, beta=0, at iter=" + str(iter) + " with residual=" + str(sqrt(inner(r,r))), 0)
#            return (x,iter,sqrt(inner(r,r)))
            return x 
        rr    = rrn
        p     = r+beta*(p-w*Ap)
        iter += 1
        #debug("r=",sqrt(inner(r,r)))
    debug("BiCGStab finished, iter: " + str(iter) + "|r|= " + str(sqrt(inner(r,r))), 0)
    if info:
        info = {}
        info["iter"] = iter
        print "Her jeg her" 
        return x, info
    return x


def precondBiCGStab(B, A, x, b, tolerance=1.0E-05, relativeconv=False, maxit=200):
    """
    precondBiCGStab(B,A,x,b): Solve Ax = b with the preconditioned biconjugate
    gradient stabilized method.

    @param B: Preconditioner supporting the __mul__ operator for operating on
    x.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """


    r = b - A*x

    debug("b0=%e"%sqrt(inner(b,b)),1)
    p = r.copy()
    rs= r.copy()
    rr = inner(r,rs)



    if relativeconv:
        tolerance *= sqrt(inner(B*r,r))
    iter = 0
    # Alloc w
    debug("r0=%e"%sqrt(inner(r,r)),1)
    while sqrt(inner(r,r)) > tolerance and iter <= maxit:
        #debug("iter, sqrt(inner(r,r))", iter, sqrt(inner(r,r)))
        ph    = B*p
        Ap    = A*ph
        alpha = rr/inner(rs,Ap)
        s     = r-alpha*Ap
        sh    = B*s
        As    = A*sh
        w     = inner(As,s)/inner(As,As)
        x    += alpha*ph+w*sh
        r     = s -w*As
        rrn   = inner(r,rs)
        beta  = (rrn/rr)*(alpha/w)
        if beta==0.0:
            debug("BiCGStab breakdown, beta=0, at iter=" + str(iter)+" with residual=" + str(sqrt(inner(r,r))), 0)
#            return (x,iter,sqrt(inner(r,r)))
            return x
        debug("|r|_%d = %e " %(iter,sqrt(inner(r,r))), 1)
        rr    = rrn
        p     = r+beta*(p-w*Ap)
        iter += 1

    debug("precondBiCGStab finished, iter: %d, ||e||=%e" % (iter,sqrt(inner(r,r))),1)
    return (x,iter,sqrt(inner(r,r)))


def precondLBiCGStab(B, A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    precondBiCGStab(B,A,x,b): Solve Ax = b with the preconditioned biconjugate
    gradient stabilized method.

    @param B: Preconditioner supporting the __mul__ operator for operating on
    x.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """

    r = b - A*x
    p = 1.0*r
    rs= 1.0*r
    rr = inner(r,rs)

    if relativeconv:
        tolerance *= sqrt(inner(b,b))
    iter = 0
    while sqrt(inner(r,r)) > tolerance:
        Ap    = A*p
        BAp   = B*Ap
        alpha = rr/inner(rs,BAp)
        s     = r-alpha*BAp
        As    = A*s
        BAs   = B*As
        w     = inner(BAs,s)/inner(BAs,BAs)
        x    += alpha*p+w*s
        r     = s -w*BAs
        rrn   = inner(r,rs)
        beta  = (rrn/rr)*(alpha/w)
        rr    = rrn
        p     = r+beta*(p-w*BAp)
        iter += 1
    debug("precondBiCGStab finished, iter: %d, ||e||=%e" % (iter,sqrt(inner(r,r))))
    return (x,iter)


"""Conjugate Gradients Method for the normal equations"""

def CGN_AA(A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    CGN_AA(B,A,x,b): Solve Ax = b with the conjugate
    gradient method applied to the normal equation.


    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """

    """
    Allocate some memory
    """
    r_a   = 1.0*x
    p_aa  = 1.0*x


    r = b - A*x
    bnrm2=sqrt(inner(b,b))
    A.prod(r,r_a,True)

    rho=inner(r_a,r_a)
    rho1=rho
    p = 1.0*r_a

    # Used to compute an estimate of the condition number
    alpha_v=[]
    beta_v=[]
    iter=0
    error   = sqrt(rho)/bnrm2
    debug("error0="+str(error), 0)

    while error>=tolerance:
        p_a     = A*p
        A.prod(p_a,p_aa)
        alpha   = rho/inner(p,p_aa)
        x       = x + alpha*p
        r       = b-A*x
        A.prod(r,r_a,True)
        rho     = inner(r_a,r_a)
        beta    = rho/rho1
        error   = sqrt(rho)/bnrm2
        debug("error = " + str(error), 0)

        alpha_v.append(alpha)
        beta_v.append(beta)

        rho1    = rho
        p       = r_a+beta*p
        iter   += 1


    debug("CGN_ABBA finished, iter: %d, ||e||=%e" % (iter,error))
    return (x,iter,alpha_v,beta_v)

def CGN_ABBA(B, A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    CGN_ABBA(B,A,x,b): Solve Ax = b with the preconditioned conjugate
    gradient method applied to the normal equation.

    @param B: Preconditioner supporting the __mul__ operator for operating on
    x.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """

    """
    Allocate some memory
    """
    r_bb   = 1.0*x
    r_abb  = 1.0*x
    p_bba  = 1.0*x
    p_abba = 1.0*x

    r = b - A*x

    bnrm2 = sqrt(inner(b,b))
    r_b = B*r
    r_bb = B*r_b #Should be transposed

    A.data[0][0].prod(r_bb[0],r_abb[0])
    A.prod(r_bb,r_abb,True)

    rho=inner(r_abb,r_abb)
    rho1=rho
    p=r_abb.copy()

    # Used to compute an estimate of the condition number
    alpha_v=[]
    beta_v=[]
    iter=0
    error  = sqrt(rho)/bnrm2

    while error>=tolerance:
        p_a     = A*p
        p_ba    = B*p_a;
        p_bba   = B*p_ba
        A.prod(p_bba,p_abba)
        alpha   = rho/inner(p,p_abba)
        x       = x + alpha*p
        r       = b-A*x
        r_b     = B*r
        r_bb    = B*r_b
        A.prod(r_bb,r_abb,True)
        rho     = inner(r_abb,r_abb)
        beta    = rho/rho1
        error   = sqrt(rho)/bnrm2
        debug("i = ", error)

        alpha_v.append(alpha)
        beta_v.append(beta)

        rho1    = rho
        p       = r_abb+beta*p
        iter   += 1


    debug("CGN_ABBA finished, iter: %d, ||e||=%e" % (iter,error))
    return (x,iter,alpha_v,beta_v)




def CGN_BABA(B, A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    CGN_BABA(B,A,x,b): Solve Ax = b with the preconditioned conjugate
    gradient method applied to the normal equation.

    @param B: Preconditioner supporting the __mul__ operator for operating on
    x.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """

    # Is this correct?? Should we not transpose the preconditoner somewhere??

    """
    Allocate some memory
    """
    r_ab   = 1.0*x
    p_aba  = 1.0*x

    r      = b - A*x
    bnrm2  = sqrt(inner(b,b))
    r_b    = B*r
    A.prod(r_b,r_ab,True)
    r_bab  = B*r_ab

    rho     = inner(r_bab,r_ab)
    rho1    = rho
    p       = r_bab.copy()

    # Used to compute an estimate of the condition number
    alpha_v=[]
    beta_v=[]
    iter=0
    error   = sqrt(rho)/bnrm2

    while error>=tolerance:
        p_a     = A*p
        p_ba    = B*p_a;
        A.prod(p_ba,p_aba,True)
        p_baba  = B*p_aba
        alpha   = rho/inner(p,p_aba)
        x       = x + alpha*p
        r       = b-A*x
        r_b     = B*r
        A.prod(r_b,r_ab,True)
        r_bab   = A*r_ab
        rho     = inner(r_bab,r_ab)
        beta    = rho/rho1
        error   = sqrt(rho)/bnrm2

        alpha_v.append(alpha)
        beta_v.append(beta)

        rho1    = rho
        p       = r_bab+beta*p
        iter   += 1


    debug("CGN_BABA finished, iter: %d, ||e||=%e" % (iter,error))
    return (x,iter,alpha_v,beta_v)



def precondRconjgrad(B, A, x, b, tolerance=1.0E-05, relativeconv=False):
    """
    conjgrad(B,A,x,b): Solve Ax = b with the right preconditioned conjugate gradient method.

    @param B: Preconditioner supporting the __mul__ operator for operating on
    x.

    @param A: Matrix or discrete operator. Must support A*x to operate on x.

    @param x: Anything that A can operate on, as well as support inner product
    and multiplication/addition/subtraction with scalar and something of the
    same type

    @param b: Right-hand side, probably the same type as x.

    @param tolerance: Convergence criterion
    @type tolerance: float

    @param relativeconv: Control relative vs. absolute convergence criterion.
    That is, if relativeconv is True, ||r||_2/||r_init||_2 is used as
    convergence criterion, else just ||r||_2 is used.  Actually ||r||_2^2 is
    used, since that save a few cycles.

    @type relativeconv: bool

    @return:  the solution x.

    """

    r0 = b - A*x
    r = 1.0*r0
    z = B*r
    q = z.copy()
    rho = rho0 = inner(z, r)
    if relativeconv:
        tolerance *= sqrt(inner(B*b,b))
    alpha_v=[]
    beta_v=[]
    iter=0

    while sqrt(fabs(rho0)) > tolerance:
        Aq=A*q
        alpha = rho0/inner(Aq,q)

        x += alpha*q
        r -= alpha*Aq
        z  = B*r
        rho = inner(z, r)
        beta = rho/rho0
        q = z+beta*q
        rho0 = rho
        alpha_v.append(alpha)
        beta_v.append(beta)
        iter += 1
    return (x,iter,alpha_v,beta_v)


def Richardson(A, x, b, tau=1, tolerance=1.0E-05, relativeconv=False, maxiter=1000, info=False):

    r = b - A*x
    rho = rho0 = inner(r, r)
    if relativeconv:
        tolerance *= sqrt(inner(b,b))
    print "tolerance ", tolerance
    iter = 0
    while sqrt(rho) > tolerance and iter <= maxiter:
        x += tau*r
        r = b - A*x
        rho = inner(r,r)
        print "iteration ", iter, " rho= ",  sqrt(rho)
        iter += 1
    return x 


def precRichardson(B, A, x, b, tau=1, tolerance=1.0E-05, relativeconv=False, maxiter=1000, info=False):

    r = b - A*x
    s = B*r
    rho = rho0 = inner(r, r)
    if relativeconv:
        tolerance *= sqrt(inner(b,b))
    print "tolerance ", tolerance
    iter = 0
    while sqrt(rho) > tolerance and iter <= maxiter:
        x += tau*s
        r = b - A*x
        s = B*r 
        rho = inner(r,r)
        print "iteration ", iter, " rho= ",  sqrt(rho)
        iter += 1
    return x 






