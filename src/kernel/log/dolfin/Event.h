// Copyright (C) 2004-2005 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2004-01-03
// Last changed: 2005

#ifndef __EVENT_H
#define __EVENT_H

#include <string>

namespace dolfin
{

  /// A event is a string message which is displayed
  /// only a limited number of times.
  ///
  /// Example of usage:
  ///
  ///   Event event("System is stiff, damping is needed.");
  ///   while ()
  ///   {
  ///     ...
  ///     if ( ... )
  ///     {
  ///       event();
  ///       ...
  ///     }
  ///   }

  class Event
  {
  public:

    /// Constructor
    Event(const std::string message, unsigned int maxcount = 1);

    /// Destructor
    ~Event();

    /// Display message
    void operator() ();
    
    /// Display count
    unsigned int count() const;

    /// Maximum display count
    unsigned int maxcount() const;    

  private:

    std::string message;
    unsigned int _maxcount;
    unsigned int _count;

  };

}

#endif
