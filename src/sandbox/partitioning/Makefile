PARMETIS = /usr/src/ParMetis-3.1
#SCOTCH = $(HOME)/shared/apps/linux/scotch_4.0/bin
SCOTCH = $(HOME)/shared/apps/linux/scotch_5.0/bin
#SCOTCH = /usr/src/scotch_4.0

CFLAGS = -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_VERSION=4.0 -DHAVE_SCOTCH_H=1 -DDEBUG -g -O2 -Wall -ansi -DNDEBUG -pedantic -Wno-long-long -std=c++98 -I../../../local/include -I/usr/include/libxml2 -I$(PARMETIS) -I$(PARMETIS)/METISLib -I$(SCOTCH)
#CFLAGS  = `pkg-config --cflags dolfin` -I$(PARMETIS) -I$(PARMETIS)/METISLib # MPI gives warnings that are being treated as errors
LIBS    = `pkg-config --libs dolfin` -L$(PARMETIS) -L$(SCOTCH) -lparmetis -lmetis -lscotch -lscotcherrcom -lcommon
CXX     = `pkg-config --variable=compiler dolfin`

DEST    = benchmark
OBJECTS = main.o

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	 $(CXX) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CFLAGS) -c $<

