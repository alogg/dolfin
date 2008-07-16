// Copyright (C) 2008 Dag Lindbo
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-07-06
// Last changed: 2008-07-06

#ifndef __MTL4_VECTOR_H
#define __MTL4_VECTOR_H

#ifdef HAS_MTL4

#include <dolfin/log/LogStream.h>
#include <dolfin/common/Variable.h>
#include "mtl4.h"
#include "GenericVector.h"

namespace dolfin
{

  class MTL4Vector: public GenericVector, public Variable
  {
  public:

    /// Create empty vector
    MTL4Vector();

    /// Create vector of size N
    explicit MTL4Vector(uint N);

    /// Copy constructor
    explicit MTL4Vector(const MTL4Vector& x);

    /// Create vector view from given MTL4_vector pointer
    explicit MTL4Vector(MTL4_vector* vector);

    /// Destructor
    virtual ~MTL4Vector();

    //--- Implementation of the GenericTensor interface ---

    /// Return copy of tensor
    virtual MTL4Vector* copy() const;

    /// Set all entries to zero and keep any sparse structure
    virtual void zero();

    /// Finalize assembly of tensor
    virtual void apply();

    /// Display vector
    virtual void disp(uint precision=2) const;

    //--- Implementation of the GenericVector interface ---

    /// Initialize vector of size N
    virtual void init(uint N);

    /// Return size of vector
    virtual uint size() const;

    /// Get block of values
    virtual void get(real* block, uint m, const uint* rows) const;

    /// Set block of values
    virtual void set(const real* block, uint m, const uint* rows);

    /// Add block of values
    virtual void add(const real* block, uint m, const uint* rows);

    /// Get all values
    virtual void get(real* values) const;

    /// Set all values
    virtual void set(real* values);

    /// Add all values to each entry
    virtual void add(real* values);

    /// Add multiple of given vector (AXPY operation)
    virtual void axpy(real a, const GenericVector& x);

    /// Return inner product with given vector
    virtual real inner(const GenericVector& vector) const;

    /// Return norm of vector
    virtual real norm(VectorNormType type = l2) const;

    /// Return minimum value of vector
    virtual real min() const;

    /// Return maximum value of vector
    virtual real max() const;

    /// Multiply vector by given number
    virtual const MTL4Vector& operator*= (real a);

    /// Divide vector by given number
    virtual const MTL4Vector& operator/= (real a)
    { *this *= 1.0 / a; return *this; }

    /// Add given vector
    virtual const MTL4Vector& operator+= (const GenericVector& x);

    /// Subtract given vector
    virtual const MTL4Vector& operator-= (const GenericVector& x);

    /// Assignment operator
    virtual const MTL4Vector& operator= (const GenericVector& x);

    /// Assignment operator
    virtual const MTL4Vector& operator= (real a);

    //--- Special functions ---
    virtual LinearAlgebraFactory& factory() const;

    //--- Special MTL4 functions ---

    /// Return MTL4_vector reference
    MTL4_vector& vec() const;

    /// Assignment operator
    const MTL4Vector& operator= (const MTL4Vector& x);

    friend class MTL4Matrix;

  private:

    // MTL4_vector pointer
    MTL4_vector* x;

    // True if we don't own the vector x points to
    bool is_view;

  };  

  LogStream& operator<< (LogStream& stream, const MTL4Vector& A);

}

#endif 
#endif 
