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

DEBUGLIBS1D = $(CLASSDIR)/graphics/1d/libgraphics.a $(CLASSDIR)/memdebug/d/libmemdebug.a -lgthread-2.0
OPTLIBS1D = $(CLASSDIR)/graphics/1o/libgraphics.a $(CLASSDIR)/memdebug/o/libmemdebug.a
DEBUGLIBS2D = $(CLASSDIR)/graphics/2d/libgraphics.a $(CLASSDIR)/memdebug/d/libmemdebug.a -lgthread-2.0
OPTLIBS2D = $(CLASSDIR)/graphics/2o/libgraphics.a $(CLASSDIR)/memdebug/o/libmemdebug.a
DEBUGLIBS3D = $(CLASSDIR)/graphics/3d/libgraphics.a $(CLASSDIR)/memdebug/d/libmemdebug.a -lgthread-2.0
OPTLIBS3D = $(CLASSDIR)/graphics/3o/libgraphics.a $(CLASSDIR)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(CLASSDIR)/memdebug -I$(CLASSDIR)/graphics `pkg-config --cflags gtkglext-x11-1.0`

#for debugging:
COMPILE.C++ = $(C++) $(DEBUGFLAGS) -DTIMING -DMEM_DEBUG -DDEBUG -DUSE_GTK
COMPILE.f = $(FC) $(DEBUGFFLAGS)

#for optimizing:
#COMPILE.C++ = $(C++) $(OPTFLAGS) -DTIMING -DUSE_GTK
#COMPILE.f = $(FC) $(OPTFLAGS)

%.o: %.C 
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<

all : lp_standard_form_objective kuhn_tucker kuhn_tucker2 penalty \
  interiorPoint sqp
lp_standard_form_objective : FORCE
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD -DSPACEDIM=3 $(INCLUDE_OPTION) -c -o lp_standard_form_objective.o lp_standard_form_objective.C
	$(C++) -o $@ lp_standard_form_objective.o $(DEBUGLIBS3D) $(GLIBS)
kuhn_tucker : FORCE
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o kuhn_tucker.o kuhn_tucker.C
	$(C++) -o $@ kuhn_tucker.o $(DEBUGLIBS2D) $(GLIBS)
kuhn_tucker2 : FORCE
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o kuhn_tucker2.o kuhn_tucker2.C
	$(C++) -o $@ kuhn_tucker2.o $(DEBUGLIBS2D) $(GLIBS)
penalty : FORCE
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o penalty.o penalty.C
	$(C++) -o $@ penalty.o $(DEBUGLIBS2D) $(GLIBS)
interiorPoint : FORCE
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o interiorPoint.o interiorPoint.C
	$(C++) -o $@ interiorPoint.o $(DEBUGLIBS2D) $(GLIBS)
sqp : FORCE
	$(COMPILE.C++) -DUSE_GTK -DNO_THREAD -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o sqp.o sqp.C
	$(C++) -o $@ sqp.o $(DEBUGLIBS2D) $(GLIBS)
#	$(MAKE) depends.gnu
clean : 
	-rm *.o core*
	-rm lp_standard_form_objective kuhn_tucker kuhn_tucker2 penalty \
	  interiorPoint sqp
depends.gnu : FORCE
	-$(RM) depends.gnu
	-cat *.d > depends.gnu
FORCE :

#include depends.gnu
