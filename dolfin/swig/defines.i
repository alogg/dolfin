/* -*- C -*- */
// Copyright (C) 2009 Johan Hake
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN.  If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2009-09-03
// Last changed: 2010-01-15

// ===========================================================================
// SWIG directives for mapping defines to Python
// ===========================================================================

%inline %{
bool has_mpi()
{
#ifdef HAS_MPI
  return true;
#else
  return false;
#endif
}

bool has_slepc()
{
#ifdef HAS_SLEPC
  return true;
#else
  return false;
#endif
}

bool has_trilinos()
{
#ifdef HAS_TRILINOS
  return true;
#else
  return false;
#endif
}

bool has_scotch()
{
#ifdef HAS_SCOTCH
  return true;
#else
  return false;
#endif
}

bool has_cgal()
{
#ifdef HAS_CGAL
  return true;
#else
  return false;
#endif
}

bool has_umfpack()
{
#ifdef HAS_UMFPACK
  return true;
#else
  return false;
#endif
}

bool has_cholmod()
{
#ifdef HAS_CHOLMOD
  return true;
#else
  return false;
#endif
}

bool has_parmetis()
{
#ifdef HAS_PARMETIS
  return true;
#else
  return false;
#endif
}

bool has_gmp()
{
#ifdef HAS_GMP
  return true;
#else
  return false;
#endif
}

bool has_zlib()
{
#ifdef HAS_ZLIB
  return true;
#else
  return false;
#endif
}

// ---------------------------------------------------------------------------
// Define a function that return true; if a specific la backend is supported
// ---------------------------------------------------------------------------
bool has_la_backend(std::string backend)
{
  if (backend == "uBLAS")
  {
    return true;
  }
  else if (backend == "PETSc")
  {
#ifdef HAS_PETSC
    return true;
#else
    return false;
#endif
  }
  else if (backend == "Epetra")
  {
#ifdef HAS_TRILINOS
    return true;
#else
    return false;
#endif
  }
  else if (backend == "MTL4")
  {
#ifdef HAS_MTL4
    return true;
#else
    return false;
#endif
  }
  else if (backend == "STL")
  {
    return true;
  }
  return false;
}

%}

%feature("docstring") has_linear_algebra_backend "
Returns True if a linear algebra backend is available.
";

