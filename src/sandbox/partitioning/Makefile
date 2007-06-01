PARMETIS=/usr/src/ParMetis-3.1/

#CFLAGS  = `pkg-config --cflags dolfin` -I$(PARMETIS) -I$(PARMETIS)/METISLib
CFLAGS  = `pkg-config --cflags dolfin` -I$(PARMETIS)/METISLib
#LIBS    = `pkg-config --libs dolfin` -L$(PARMETIS) -L$(SCOTCH_LIB) -lparmetis -lmetis -lmpi -llapack -lX11 -lm
LIBS    = `pkg-config --libs dolfin` -L$(SCOTCH_LIB) -lmetis -llapack -lX11 -lm
CXX     = `pkg-config --variable=compiler dolfin`

DEST    = dolfin-partition
OBJECTS = main.o

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	 $(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
