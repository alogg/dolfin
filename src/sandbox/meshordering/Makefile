CFLAGS  = `pkg-config --cflags dolfin`
LIBS    = `pkg-config --libs dolfin`
CXX     = `pkg-config --variable=compiler dolfin`

DEST    = dolfin-meshordering
#OBJECTS = main.o MeshOrdering.o
#OBJECTS = main.o Triangle.o
#OBJECTS = main.o Tetrahedron.o
OBJECTS = main.o

all: $(DEST)

install:

clean:
	-rm -f *~ *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	 $(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
