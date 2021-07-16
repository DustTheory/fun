SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
LAPACKPP_BASE = $(DIR_BASE)/lapack++
DEBUGLIBS = $(DIR_BASE)/gui/d/libgui.a \
  $(DIR_BASE)/graphics/1d/libgraphics.a \
  $(DIR_BASE)/memdebug/d/libmemdebug.a \
  -lgthread-2.0
OPTLIBS = $(DIR_BASE)/gui/o/libgui.a \
  $(DIR_BASE)/graphics/1o/libgraphics.a \
  $(DIR_BASE)/memdebug/o/libmemdebug.
INCLUDE_OPTION = -I$(SD) -I$(DIR_BASE)/gui -I$(DIR_BASE)/graphics \
  -I$(DIR_BASE)/memdebug -I$(LAPACKPP_BASE)\
  `pkg-config --cflags gtkglext-x11-1.0`
GLIBS = `pkg-config --libs gtkglext-x11-1.0`

FC = gfortran
DEBUGFFLAGS = -c -g -fimplicit-none -fno-f2c -finit-real=INF -fsecond-underscore   -Werror
FLIBS = -lgfortran

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
DEBUGFLAGS = -g
OPTFLAGS = -O3 -g

CC = gcc -MMD

#for debugging:
COMPILE.C++ = $(C++) $(DEBUGFLAGS) -DTIMING -DMEM_DEBUG -DDEBUG -DUSE_GTK -DNO_THREAD
COMPILE.f = $(FC) $(DEBUGFFLAGS)
CLASS_LIBS = $(DEBUGLIBS)
LAPACKPP_LIB = $(LAPACKPP_BASE)/d/lapack++.a
#for optimizing:
#COMPILE.C++ = $(C++) $(OPTFLAGS) -DTIMING -DUSE_GTK
#COMPILE.f = $(FC) $(OPTFLAGS)
#CLASS_LIBS = $(OPTLIBS)
#LAPACKPP_LIB = $(LAPACKPP_BASE)/o/lapack++.a

%.o: %.C
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<
all : GUIFiniteDifference GUIFiniteElement
GUIFiniteDifference : GUIFiniteDifference.o finite_difference.o
	$(C++) -o $@ GUIFiniteDifference.o finite_difference.o \
	$(CLASS_LIBS) $(LIBS) $(SLIBS) $(GLIBS) $(FLIBS)
	$(MAKE) depends.gnu

GUIFiniteElement : GUIFiniteElement.o Polynomial.o Quadrature.o
	$(C++) -o $@ GUIFiniteElement.o Polynomial.o Quadrature.o \
	$(CLASS_LIBS) $(LAPACKPP_LIB) $(SLIBS) $(GLIBS) $(FLIBS) \
	-llapack -lblas
	$(MAKE) depends.gnu

#GUIFiniteElement3 : GUIFiniteElement3.o finite_element.o
#	$(C++) -o $@ GUIFiniteElement3.o finite_element.o \
#	$(CLASS_LIBS) $(LIBS) $(SLIBS) -llapack
#	$(MAKE) depends.gnu

#GUIHierarchical : GUIHierarchical.o hierarchical.o
#	$(C++) -o $@ GUIHierarchical.o hierarchical.o \
#	$(CLASS_LIBS) $(LIBS) $(SLIBS) -llapack
#	$(MAKE) depends.gnu

clean : FORCE
	-rm *.o trash* core*
	-rm GUIFiniteDifference GUIFiniteElement
FORCE :

depends.gnu :
	-$(RM) depends.gnu
	cat *.d > depends.gnu
include depends.gnu
