// Copyright (C) 2003 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#ifndef __DIAGONAL_ITERATION_H
#define __DIAGONAL_ITERATION_H

#include <dolfin/NewArray.h>
#include <dolfin/Iteration.h>

namespace dolfin
{

  /// State-specific behavior of fixed point iteration for diagonally stiff problems.

  class DiagonalIteration : public Iteration
  {
  public:

    DiagonalIteration(Solution& u, RHS& f, FixedPointIteration& fixpoint,
		      unsigned int maxiter, real maxdiv, real maxconv, real tol);
    
    ~DiagonalIteration();

    State state() const;

    void start(TimeSlab& timeslab);
    void start(NewArray<Element*>& elements);
    void start(Element& element);

    void update(TimeSlab& timeslab);
    void update(NewArray<Element*>& elements);
    void update(Element& element);    

    void stabilize(TimeSlab& timeslab, const Residuals& r, unsigned int n);
    void stabilize(NewArray<Element*>& elements, const Residuals& r, unsigned int n);
    void stabilize(Element& element, const Residuals& r, unsigned int n);
    
    bool converged(TimeSlab& timeslab, Residuals& r, unsigned int n);
    bool converged(NewArray<Element*>& elements, Residuals& r, unsigned int n);
    bool converged(Element& element, Residuals& r, unsigned int n);

    bool diverged(TimeSlab& timeslab, Residuals& r, unsigned int n, Iteration::State& newstate);
    bool diverged(NewArray<Element*>& elements, Residuals& r, unsigned int n, Iteration::State& newstate);
    bool diverged(Element& element, Residuals& r, unsigned int n, Iteration::State& newstate);

    void report() const;

  private:

    // Compute alpha
    real computeAlpha(real rho) const;

    // Stabilization parameter
    real alpha;

  };

}

#endif
