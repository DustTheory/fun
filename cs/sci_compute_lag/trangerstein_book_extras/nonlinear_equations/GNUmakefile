SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
CLASSDIR = $(DIR_BASE)

DEBUGFLAGS = -g -fno-implicit-templates -fno-default-inline
OPTFLAGS = -O3
C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -g -MMD
#LIBS = -lf2c 
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -Wall
CC = gcc -g -MMD 
FC = gfortran 
GLIBS = `pkg-config --libs gtkglext-x11-1.0`
GSLLIBS = /usr/lib64/libgsl.so

DEBUGLIBS = $(CLASSDIR)/graphics/1d/libgraphics.a $(CLASSDIR)/memdebug/d/libmemdebug.a -lgthread-2.0
OPTLIBS = $(CLASSDIR)/graphics/1o/libgraphics.a $(CLASSDIR)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(CLASSDIR)/memdebug -I$(CLASSDIR)/graphics `pkg-config --cflags gtkglext-x11-1.0`

#for debugging:
COMPILE.C++ = $(C++) $(DEBUGFLAGS) -DTIMING -DMEM_DEBUG -DDEBUG -DUSE_GTK
COMPILE.f = $(FC) $(DEBUGFFLAGS)
LIBS = $(DEBUGLIBS)

#for optimizing:
#COMPILE.C++ = $(C++) $(OPTFLAGS) -DTIMING -DUSE_GTK
#COMPILE.f = $(FC) $(OPTFLAGS)
#LIBS = $(OPTLIBS)

%.o: %.C 
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<

all : bisection newton fdnewton secant efficiency illinois mullerTraub rationalInterpolation inverseQuadraticInterpolation parabolicMin
bisection : bisectionMain.o
	$(C++) -o $@ bisectionMain.o $(LIBS) $(GLIBS) $(GSLLIBS)
	$(MAKE) depends.gnu
newton : newtonMain.o
	$(C++) -o $@ newtonMain.o $(LIBS) $(GLIBS) $(GSLLIBS)
	$(MAKE) depends.gnu
fdnewton : finiteDifferenceNewton.o
	$(C++) -o $@ finiteDifferenceNewton.o $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu
secant : secantMain.o
	$(C++) -o $@ secantMain.o $(LIBS) $(GLIBS) $(GSLLIBS)
	$(MAKE) depends.gnu
efficiency : relativeEfficiency.o
	$(C++) -o $@ relativeEfficiency.o $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu
#regulaFalsi : regulaFalsiMain.o
#	$(C++) -o $@ regulaFalsiMain.o $(LIBS) $(GLIBS)
#	$(MAKE) depends.gnu
illinois : illinoisMain.o
	$(C++) -o $@ illinoisMain.o $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu
mullerTraub : mullerTraub.o
	$(C++) -o $@ mullerTraub.o $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu
rationalInterpolation : rationalInterpolation.o
	$(C++) -o $@ rationalInterpolation.o $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu
inverseQuadraticInterpolation : inverseQuadraticInterpolation.o
	$(C++) -o $@ inverseQuadraticInterpolation.o $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu
parabolicMin : parabolicMin.o fmin.o d1mach.o
	$(C++) -o $@ parabolicMin.o fmin.o d1mach.o $(LIBS) $(GLIBS) -lgfortran
	$(MAKE) depends.gnu
depends.gnu : FORCE
	-$(RM) depends.gnu
	cat *.d > depends.gnu
clean : 
	-rm *.o core*
	-rm bisection newton fdnewton secant efficiency illinois mullerTraub
	-rm rationalInterpolation inverseQuadraticInterpolation parabolicMin
FORCE :

#include depends.gnu
