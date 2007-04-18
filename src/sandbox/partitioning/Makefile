PARMETIS=/home/magnus/doc/uio/master/ParMetis-3.1/
SCOTCH_LIB=/home/magnus/shared/apps/linux/scotch_4.0/bin
SCOTCH=/home/magnus/shared/apps/linux/scotch_4.0/src/libscotch

CFLAGS  = `pkg-config --cflags dolfin` -I$(PARMETIS) -I$(PARMETIS)/METISLib -I$(SCOTCH) -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_VERSION="\"4.0\""
LIBS    = `pkg-config --libs dolfin` -L$(PARMETIS) -L$(SCOTCH_LIB) -lparmetis -lmetis -lmpi -llapack -lX11 -lm -lscotch -lscotcherrcom -lcommon
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
