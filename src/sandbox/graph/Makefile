PARMETIS=/home/magnus/doc/uio/master/ParMetis-3.1/

CFLAGS  = `pkg-config --cflags dolfin` -I$(PARMETIS) -I$(PARMETIS)/METISLib
LIBS    = `pkg-config --libs dolfin` -L$(PARMETIS) -lparmetis -lmetis
CXX     = `pkg-config --variable=compiler dolfin`

DEST    = dolfin-graph
OBJECTS = main.o

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	 $(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
