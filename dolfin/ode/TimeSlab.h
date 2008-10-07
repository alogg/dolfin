// Copyright (C) 2005-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-05-02
// Last changed: 2008-10-06

#ifndef __TIME_SLAB_H
#define __TIME_SLAB_H

#include <dolfin/common/types.h>

namespace dolfin
{

  class ODE;
  class Method;
  class uBLASVector;

  /// This is the base class for time slabs, the collections of
  /// degrees of freedom for the solution of an ODE between two
  /// synchronized time levels a and b.

  class TimeSlab
  {
  public:

    /// Constructor
    TimeSlab(ODE& ode);

    /// Destructor
    virtual ~TimeSlab();
    
    /// Build time slab, return end time
    virtual double build(double a, double b) = 0;

    /// Solve time slab system
    virtual bool solve() = 0;

    /// Check if current solution can be accepted
    virtual bool check(bool first) = 0;

    /// Shift time slab (prepare for next time slab)
    virtual bool shift(bool end) = 0;

    /// Prepare sample at time t
    virtual void sample(double t) = 0;

    /// Return number of components
    uint size() const;

    /// Return start time of time slab
    double starttime() const;
    
    /// Return end time of time slab
    double endtime() const;

    /// Return length of time slab
    double length() const;

    /// Sample solution value of given component at given time
    virtual double usample(uint i, double t) = 0;

    /// Sample time step size for given component at given time
    virtual double ksample(uint i, double t) = 0;

    /// Sample residual for given component at given time
    virtual double rsample(uint i, double t) = 0;

    /// Display time slab data
    virtual void disp() const = 0;

    /// Output
    friend LogStream& operator<<(LogStream& stream, const TimeSlab& timeslab);

    /// Friends
    friend class TimeSlabJacobian;
    friend class TimeSlabSolver;

  protected:

    // Write given solution vector to file
    static void write(uint N, const double* u);

    // Copy data of given size between vectors with given offsets
    static void copy(const double* x, uint xoffset, double* y, uint yoffset, uint n);

    // Copy data of given size between vectors with given offsets
    static void copy(const uBLASVector& x, uint xoffset, double* y, uint yoffset, uint n);

    // Copy data of given size between vectors with given offsets
    static void copy(const double* x, uint xoffset, uBLASVector& y, uint yoffset, uint n);

    // Copy data of given size between vectors with given offsets
    static void copy(const uBLASVector& x, uint xoffset, uBLASVector& y, uint yoffset, uint n);
    
    uint N;    // Size of system
    double _a; // Start time of time slab
    double _b; // End time of time slab
    
    ODE& ode;             // The ODE
    const Method* method; // Method, mcG(q) or mdG(q)  
    double* u0;           // Initial values
    
    bool save_final; // True if we should save the solution at final time

  };

}

#endif
