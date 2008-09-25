// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
//
// First added:  2007-04-02
// Last changed: 2008-09-25

#ifndef __FORM_H
#define __FORM_H

#include <tr1/memory>
#include "DofMapSet.h"

// Forward declaration
namespace ufc 
{
  class form; 
}

namespace dolfin
{

  // Forward declarations
  class Function;
  template<class T> class Array;
  template<class T> class MeshFunction;


  /// Base class for UFC code generated by FFC for DOLFIN with option -l

  class Form
  {
  public:

    /// Constructor
    Form() {}

    /// Destructor
    virtual ~Form();

    /// Return UFC form
    virtual const ufc::form& form() const = 0;

    /// Return array of coefficients
    virtual const Array<Function*>& coefficients() const = 0;

    /// Create degree of freedom maps 
    void updateDofMaps(Mesh& mesh);

    /// Create degree of freedom maps
    void updateDofMaps(Mesh& mesh, MeshFunction<uint>& partitions);

    /// Set degree of freedom maps (not owned)
    void setDofMaps(DofMapSet& dof_map_set);

    /// Return DofMapSet
    DofMapSet& dofMaps() const;

  private:

    // Degree of freedom maps
    std::tr1::shared_ptr<DofMapSet> dof_map_set;

  };

}

#endif
