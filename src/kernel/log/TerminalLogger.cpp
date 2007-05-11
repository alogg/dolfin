// Copyright (C) 2003-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2003-03-13
// Last changed: 2007-05-11

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <dolfin/utils.h>
#include <dolfin/log.h>
#include <dolfin/TerminalLogger.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
TerminalLogger::TerminalLogger() : GenericLogger()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
TerminalLogger::~TerminalLogger()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void TerminalLogger::info(const char* msg)
{
  indent();
  std::cout << msg << std::endl;
}
//-----------------------------------------------------------------------------
void TerminalLogger::debug(const char* msg, const char* location)
{
  indent();
  std::cout << "[Debug at " << location << "]: " << msg << std::endl;
}
//-----------------------------------------------------------------------------
void TerminalLogger::warning(const char* msg, const char* location)
{
  indent();
  std::cout << "*** Warning: " << msg << std::endl;
}
//-----------------------------------------------------------------------------
void TerminalLogger::error(const char* msg, const char* location)
{
  indent();
  std::cout << "*** Error: " << msg << " [" << location
	    << "]" << std::endl;
  exit(1);
}
//-----------------------------------------------------------------------------
void TerminalLogger::dassert(const char* msg, const char* location)
{
  indent();
  std::cout << "*** Assertion " << msg << " failed [" << location
	    << "]" <<  std::endl;
  dolfin_segfault();
}
//-----------------------------------------------------------------------------
void TerminalLogger::progress(const char* title, const char* label, real p)
{
  int N = DOLFIN_TERM_WIDTH - 15;
  int n = (int) (p*((real) N));
  
  // Print the title
  indent();
  printf("| %s", title);
  for (int i = 0; i < (N-length(title)-1); i++)
    printf(" ");
  printf("|\n");
  
  // Print the progress bar
  indent();
  printf("|");
  for (int i = 0; i < n; i++)
    printf("=");
  if ( n > 0 && n < N ) {
    printf("|");
    n++;
  }
  for (int i = n; i < N; i++)
    printf("-");
  printf("| %.1f%%\n", 100.0*p);
}
//-----------------------------------------------------------------------------
void TerminalLogger::indent()
{
  // Indent output to indicate the level
  for (int i = 0; i < level; i++)
    std::cout << "  ";
}
//-----------------------------------------------------------------------------
