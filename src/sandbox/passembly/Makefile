#PARMETIS = /usr/local/garth/software/ParMetis-3.1
PARMETIS=/home/magnus/doc/uio/master/ParMetis-3.1/
PETSC=/home/magnus/shared/apps/linux/petsc-2.3.2-p8/include
PETSC_LIB=/home/magnus/shared/apps/linux/petsc-2.3.2-p8/lib

CFLAGS  = `pkg-config dolfin --cflags` -I$(PARMETIS) -I$(PARMETIS)/METISLib -I$(PETSC) -DMPICH_IGNORE_CXX_SEEK
LIBS    = `pkg-config dolfin --libs` -L$(PARMETIS) -L$(PETSC_LIB) -lparmetis -lmetis -lpetc
CXX     = g++
LINK    = libtool --silent --mode=link --tag CXX g++

DEST    = dolfin-parallel-test
OBJECTS = main.o

all: $(DEST)

install:

uninstall:

clean:
	-rm -f *.o core *.core $(OBJECTS) $(DEST)

$(DEST): $(OBJECTS)
	 $(LINK) -o $@ $(OBJECTS) $(CFLAGS) $(LIBS) 

.cpp.o:
	$(CXX) $(CFLAGS) -c $<
