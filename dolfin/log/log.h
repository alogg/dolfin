// Copyright (C) 2003-2008 Anders Logg and Jim Tilander.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Ola Skavhaug, 2007, 2009.
//
// First added:  2003-03-13
// Last changed: 2009-07-02

#ifndef __LOG_H
#define __LOG_H

#include <string>
#include <cassert>
#include <dolfin/common/types.h>

namespace dolfin
{

  class Variable;
  class Parameters;

  /// The DOLFIN log system provides the following set of functions for
  /// uniform handling of log messages, warnings and errors. In addition,
  /// macros are provided for debug messages and assertions.
  ///
  /// Only messages with a debug level higher than or equal to the current
  /// log level are printed (the default being zero). Logging may also be
  /// turned off by calling log(false).

  /// Print message
  void info(std::string msg, ...);

  /// Print message at given debug level
  void info(int debug_level, std::string msg, ...);

  /// Print variable (using output of str() method)
  void info(const Variable& variable);

  /// Print variable (using output of str() method)
  void info(const Parameters& parameters);

  /// Print message to stream
  void info_stream(std::ostream& out, std::string msg);

  /// Print underlined message
  void info_underline(std:: string msg, ...);

  /// Print warning
  void warning(std::string msg, ...);

  /// Print error message and throw an exception
  void error(std::string msg, ...);

  /// Begin task (increase indentation level)
  void begin(std::string msg, ...);

  /// Begin task (increase indentation level)
  void begin(int debug_level, std::string msg, ...);

  /// End task (decrease indentation level)
  void end();
  
  /// Turn logging on or off
  void logging(bool active=true);

  /// Set log level
  void set_log_level(uint level);

  /// Get log level
  uint get_log_level();

  /// Indent string
  std::string indent(std::string s);

  /// Print summary of timings and tasks, optionally clearing stored timings
  void summary(bool reset=false);

  /// Return timing (average) for given task, optionally clearing timing for task
  double timing(std::string task, bool reset=false);

  // Helper function for dolfin_debug macro
  void __debug(std::string file, unsigned long line, std::string function, std::string format, ...);

}

// Debug macros (with varying number of arguments)
#define dolfin_debug(msg)              do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg); } while (false)
#define dolfin_debug1(msg, a0)         do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0); } while (false)
#define dolfin_debug2(msg, a0, a1)     do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1); } while (false)
#define dolfin_debug3(msg, a0, a1, a2) do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, msg, a0, a1, a2); } while (false)

#define dolfin_assert(check) assert(check)

// Not implemented error, reporting function name and line number
#define dolfin_not_implemented() do { dolfin::__debug(__FILE__, __LINE__, __FUNCTION__, "Sorry, this function has not been implemented."); error("Not implemented"); } while (false)

#endif
