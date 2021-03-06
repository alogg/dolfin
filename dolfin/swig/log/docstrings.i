// Auto generated SWIG file for Python interface of DOLFIN
//
// Copyright (C) 2012 Kristian B. Oelgaard
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//

// Autogenerated docstrings file, extracted from the DOLFIN source C++ files.

// Documentation extracted from: (module=log, header=log.h)
%feature("docstring")  dolfin::info "
**Overloaded versions**

* info\ (msg, ...)

  The DOLFIN log system provides the following set of functions for
  uniform handling of log messages, warnings and errors. In addition,
  macros are provided for debug messages and dolfin_assertions.
  
  Only messages with a debug level higher than or equal to the current
  log level are printed (the default being zero). Logging may also be
  turned off by calling set_log_active(false).
  Print message

* info\ (parameters, verbose=false)

  Print parameter (using output of str() method)

* info\ (variable, verbose=false)

  Print variable (using output of str() method)
";

%feature("docstring")  dolfin::info_stream "
Print message to stream
";

%feature("docstring")  dolfin::info_underline "
Print underlined message
";

%feature("docstring")  dolfin::warning "
Print warning
";

%feature("docstring")  dolfin::error "
Print error message and throw an exception.
Note to developers: this function should not be used internally
in DOLFIN. Use the more informative dolfin_error instead.
";

%feature("docstring")  dolfin::dolfin_error "
Print error message. Prefer this to the above generic error message.

*Arguments*
    location (str)
        Name of the file from which the error message was generated.
    task (str)
        Name of the task that failed.
        Note that this string should begin with lowercase.
        Note that this string should not be punctuated.
    reason (str)
        A format string explaining the reason for the failure.
        Note that this string should begin with uppercase.
        Note that this string should not be punctuated.
        Note that this string may contain printf style formatting.
    ... (primitive types like int, uint, double, bool)
        Optional arguments for the format string.

Developers should read the file dolfin/log/README in the DOLFIN
source tree for further notes about the use of this function.
";

%feature("docstring")  dolfin::deprecation "
Issue deprecation warning for removed feature

*Arguments*
    feature (str)
       Name of the feature that has been removed.
    version (str)
       Version number of the release in which the feature was removed.
    message (str)
       A format string explaining the deprecation.
";

%feature("docstring")  dolfin::log "
Print message at given debug level
";

%feature("docstring")  dolfin::begin "
**Overloaded versions**

* begin\ (msg, ...)

  Begin task (increase indentation level)

* begin\ (debug_level, msg, ...)

  Begin task (increase indentation level)
";

%feature("docstring")  dolfin::end "
End task (decrease indentation level)
";

%feature("docstring")  dolfin::set_log_active "
Turn logging on or off
";

%feature("docstring")  dolfin::set_log_level "
Set log level
";

%feature("docstring")  dolfin::set_output_stream "
Set output stream
";

%feature("docstring")  dolfin::get_log_level "
Get log level
";

%feature("docstring")  dolfin::monitor_memory_usage "
Monitor memory usage. Call this function at the start of a
program to continuously monitor the memory usage of the process.
";

%feature("docstring")  dolfin::not_working_in_parallel "
Report that functionality has not (yet) been implemented to work
in parallel
";

// Documentation extracted from: (module=log, header=Event.h)
%feature("docstring")  dolfin::Event "
A event is a string message which is displayed
only a limited number of times.

*Example*
    .. code-block:: python
    
        >>> event = dolfin.Event(\"System is stiff, damping is needed.\", 3)
        >>> for i in range(10):
        ...     if i > 7:
        ...         print i
        ...         event()
";

%feature("docstring")  dolfin::Event::Event "
Constructor
";

%feature("docstring")  dolfin::Event::operator "
Display message
";

%feature("docstring")  dolfin::Event::count "
Display count
";

%feature("docstring")  dolfin::Event::maxcount "
Maximum display count
";

// Documentation extracted from: (module=log, header=LogStream.h)
%feature("docstring")  dolfin::LogStream "
This class provides functionality similar to standard C++
streams (std::cout, std::endl) for output but working through
the DOLFIN log system.
";

%feature("docstring")  dolfin::LogStream::LogStream "
Create log stream of given type
";

%feature("docstring")  dolfin::LogStream::operator<< "
**Overloaded versions**

* operator<<\ (stream)

  Output for log stream

* operator<<\ (s)

  Output for string

* operator<<\ (a)

  Output for int

* operator<<\ (a)

  Output for unsigned int

* operator<<\ (a)

  Output for long int

* operator<<\ (a)

  Output for long int

* operator<<\ (a)

  Output for double

* operator<<\ (z)

  Output for std::complex<double>

* operator<<\ (variable)

  Output for variable (calling str() method)

* operator<<\ (entity)

  Output for mesh entity (not subclass of Variable for efficiency)

* operator<<\ (point)

  Output for point (not subclass of Variable for efficiency)
";

// Documentation extracted from: (module=log, header=Progress.h)
%feature("docstring")  dolfin::Progress "
This class provides a simple way to create and update progress
bars during a computation.

*Example*
    A progress bar may be used either in an iteration with a known number
    of steps:
    
    .. code-block:: python
    
        >>> n = 1000000
        >>> p = dolfin.Progress(\"Iterating...\", n)
        >>> for i in range(n):
        ...     p += 1
    
    or in an iteration with an unknown number of steps:
    
    .. code-block:: python
    
        >>> pr = dolfin.Progress(\"Iterating\")
        >>> t = 0.0
        >>> n = 1000000.0
        >>> while t < n:
        ...     t += 1.0
        ...     p += t/n
";

%feature("docstring")  dolfin::Progress::Progress "
**Overloaded versions**

* Progress\ (title, n)

  Create progress bar with a known number of steps
  
  *Arguments*
      title (str)
          The title.
      n (int)
          Number of steps.

* Progress\ (title)

  Create progress bar with an unknown number of steps
  
  *Arguments*
      title (str)
          The title.
";

%feature("docstring")  dolfin::Progress::operator= "
Set current position

*Arguments*
    p (float)
        The position.
";

%feature("docstring")  dolfin::Progress::operator++ "
Increment progress
";

// Documentation extracted from: (module=log, header=Table.h)
%feature("docstring")  dolfin::Table "
This class provides storage and pretty-printing for tables.
Example usage:

  Table table(\"Timings\");

  table(\"uBLAS\",  \"Assemble\") = 0.010;
  table(\"uBLAS\",  \"Solve\")    = 0.020;
  table(\"PETSc\",  \"Assemble\") = 0.011;
  table(\"PETSc\",  \"Solve\")    = 0.019;
  table(\"Epetra\", \"Assemble\") = 0.012;
  table(\"Epetra\", \"Solve\")    = 0.018;

  info(table);
";

%feature("docstring")  dolfin::Table::Table "
Create empty table
";

%feature("docstring")  dolfin::Table::operator "
Return table entry
";

%feature("docstring")  dolfin::Table::set "
**Overloaded versions**

* set\ (row, col, value)

  Set value of table entry

* set\ (row, col, value)

  Set value of table entry

* set\ (row, col, value)

  Set value of table entry

* set\ (row, col, value)

  Set value of table entry
";

%feature("docstring")  dolfin::Table::get "
Get value of table entry
";

%feature("docstring")  dolfin::Table::get_value "
Get value of table entry
";

%feature("docstring")  dolfin::Table::title "
Return table title
";

%feature("docstring")  dolfin::Table::operator= "
Assignment operator
";

%feature("docstring")  dolfin::Table::str "
Return informal string representation (pretty-print)
";

%feature("docstring")  dolfin::Table::str_latex "
Return informal string representation for LaTeX
";

%feature("docstring")  dolfin::TableEntry "
This class represents an entry in a Table
";

%feature("docstring")  dolfin::TableEntry::TableEntry "
Create table entry
";

%feature("docstring")  dolfin::TableEntry::operator= "
**Overloaded versions**

* operator=\ (value)

  Assign value to table entry

* operator=\ (value)

  Assign value to table entry

* operator=\ (value)

  Assign value to table entry

* operator=\ (value)

  Assign value to table entry
";

%feature("docstring")  dolfin::TableEntry::operator std::string "
Cast to entry value
";

// Documentation extracted from: (module=log, header=LogLevel.h)
