import libxml2, sys

try:
    rngfile = sys.argv[1]
    xmlfile = sys.argv[2]
except:
    print 'usage: %s RNG-file XML-file' %sys.argv[0]
    sys.exit(1)

# Read XML file:
reader = libxml2.newTextReaderFilename(xmlfile)

# Read schema file:
schema = open(rngfile).read()
rngp = libxml2.relaxNGNewMemParserCtxt(schema, len(schema))
rngs = rngp.relaxNGParse()

# Parse file:
reader.RelaxNGSetSchema(rngs)
ret = 1
while ret == 1:
    ret = reader.Read()
    if ret == -1:
        print 'Error in file %s, exiting' %xmlfile
        sys.exit(ret)

if reader.IsValid() == 1:
    print '%s validates' %xmlfile
else:
    print '%s fails to validate' %xmlfile
