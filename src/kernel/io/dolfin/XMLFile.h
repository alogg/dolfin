// Copyright (C) 2003-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2003-07-15
// Last changed: 2006-02-27

#ifndef __XML_FILE_H
#define __XML_FILE_H

#include <libxml/parser.h>

#include <dolfin/constants.h>
#include <dolfin/Vector.h>
#include <dolfin/Matrix.h>
#include <dolfin/GenericFile.h>

namespace dolfin
{
  
  class Mesh;
  class ParameterList;
  class BLASFormData;

  class XMLObject;
  
  class XMLFile : public GenericFile
  {
  public:
    
    XMLFile(const std::string filename);
    ~XMLFile();
    
    // Input

    void operator>> (Vector& x);
    void operator>> (Matrix& A);
    void operator>> (Mesh& mesh);
    void operator>> (Function& f);
    void operator>> (FiniteElementSpec& spec);
    void operator>> (ParameterList& parameters);
    void operator>> (BLASFormData& blas);
    
    void parse(Function& f, FiniteElement& element);
    
    // Output
    
    void operator<< (Vector& x);
    void operator<< (Matrix& A);
    void operator<< (Mesh& mesh);
    void operator<< (Function& f);
    void operator<< (FiniteElementSpec& spec);
    void operator<< (ParameterList& parameters);
    
    // Friends
    
    friend void sax_start_element (void *ctx, const xmlChar *name, const xmlChar **attrs);
    friend void sax_end_element   (void *ctx, const xmlChar *name);
    
  private:
    
    void parseFile();
    void parseSAX();

    FILE* openFile();
    void  closeFile(FILE* fp);
    
    // Implementation for specific class (output)
    XMLObject* xmlObject;

    // True if header is written (need to close)
    bool header_written;

    // Most recent position in file
    long mark;

  };
  
  // Callback functions for the SAX interface
  
  void sax_start_document (void *ctx);
  void sax_end_document   (void *ctx);
  void sax_start_element  (void *ctx, const xmlChar *name, const xmlChar **attrs);
  void sax_end_element    (void *ctx, const xmlChar *name);

  void sax_warning     (void *ctx, const char *msg, ...);
  void sax_error       (void *ctx, const char *msg, ...);
  void sax_fatal_error (void *ctx, const char *msg, ...);
  
}

#endif
