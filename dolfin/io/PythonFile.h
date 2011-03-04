// Copyright (C) 2003-2005 Johan Hoffman and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-07-15
// Last changed: 2005

#ifndef __PYTHON_FILE_H
#define __PYTHON_FILE_H

#include <string>
#include <tr1/tuple>
#include <dolfin/common/types.h>
#include <dolfin/common/real.h>
#include "GenericFile.h"
#include <boost/ref.hpp>

namespace dolfin
{
  typedef boost::reference_wrapper<Array<real> > RealArrayRef;

  class Sample;

  // Represents input/output of data in a format readable by Python
  // (Numeric). The data is written to several files (the variable
  // name is appended to the base file name) to enable incremental
  // output in an efficient way.

  class PythonFile : public GenericFile
  {
  public:

    PythonFile(const std::string filename);
    virtual ~PythonFile();

    // Input

    // Output

    void operator<< (const Sample& sample);
    void operator<< (const std::pair<real, RealArrayRef>);

    std::string filename_t, filename_u, filename_k, filename_r;

  };

}

#endif
