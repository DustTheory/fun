SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
CLASSDIR = $(DIR_BASE)

DEBUGFLAGS = -g -fno-implicit-templates -fno-default-inline
OPTFLAGS = -O3
C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -g -MMD
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -Wall
CC = gcc -g -MMD
FC = gfortran
GLIBS = `pkg-config --libs gtkglext-x11-1.0`
GSLLIBS = /usr/lib64/libgsl.so

DEBUGLIBS = $(CLASSDIR)/graphics/1d/libgraphics.a $(CLASSDIR)/memdebug/d/libmemdebug.a -lgthread-2.0
OPTLIBS = $(CLASSDIR)/graphics/1o/libgraphics.a $(CLASSDIR)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(CLASSDIR)/memdebug -I$(CLASSDIR)/graphics -I$(CLASSDIR)/lapack++ `pkg-config --cflags gtkglext-x11-1.0`

COMPILE.C++ = $(C++) $(DEBUGFLAGS) -DTIMING -DMEM_DEBUG -DDEBUG -DUSE_GTK
COMPILE.f = $(FC) $(FFLAGS)
LIBS = $(DEBUGLIBS)
#COMPILE.C++ = $(C++) $(OPTFLAGS) -DTIMING -DUSE_GTK
#COMPILE.f = $(FC) $(OPTFLAGS)
#LIBS = $(OPTLIBS)

%.o: %.C
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<

all : ratesConvergence mintest lennard_jones LennardJones LennardJonesSA
ratesConvergence : ratesConvergence.o
	$(LINK.C) -o ratesConvergence ratesConvergence.o $(LIBS) $(GLIBS)
mintest : testNonlinearSystemSolver.o NonlinearSystemSolver.o
	$(LINK.C) -o mintest testNonlinearSystemSolver.o \
	  NonlinearSystemSolver.o $(CLASSDIR)/lapack++/d/lapack++.a \
	  $(CLASSDIR)/minpack/libminpack.a $(LIBS) -llapack -lblas
lennard_jones : lennard_jones.o $(DIR_BASE)/minpack/libminpack.a
	$(COMPILE.f) -o $@ lennard_jones.o $(DIR_BASE)/minpack/libminpack.a
LennardJones : LennardJones.o
	$(LINK.C) -o LennardJones LennardJones.o -lgsl $(LIBS)
LennardJonesSA : LennardJonesSA.o
	$(LINK.C) -o LennardJonesSA LennardJonesSA.o -lgsl $(LIBS)
clean : 
	-rm *.o core*
	-rm ratesConvergence mintest lennard_jones LennardJones 
	-rm LennardJonesSA
