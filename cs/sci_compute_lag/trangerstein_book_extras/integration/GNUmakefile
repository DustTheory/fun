SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
DEBUGLIBS = $(DIR_BASE)/gui/d/libgui.a \
  $(DIR_BASE)/graphics/1d/libgraphics.a \
  $(DIR_BASE)/memdebug/d/libmemdebug.a \
  -lgthread-2.0
OPTLIBS = $(DIR_BASE)/gui/o/libgui.a \
  $(DIR_BASE)/graphics/1o/libgraphics.a \
  $(DIR_BASE)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(DIR_BASE)/gui -I$(DIR_BASE)/graphics \
  -I$(DIR_BASE)/memdebug `pkg-config --cflags gtkglext-x11-1.0`
GLIBS = `pkg-config --libs gtkglext-x11-1.0`

FC = gfortran
DEBUGFFLAGS = -c -g -fimplicit-none -fno-f2c -finit-real=INF -fsecond-underscore   -Wunused -Werror
FLIBS = -lgfortran

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
DEBUGFLAGS = -g
OPTFLAGS = -O3 -g

CC = gcc -MMD

#for debugging:
COMPILE.C++ = $(C++) $(DEBUGFLAGS) -DTIMING -DMEM_DEBUG -DDEBUG -DUSE_GTK -DNO_THREAD
COMPILE.f = $(FC) $(DEBUGFFLAGS)
CLASS_LIBS = $(DEBUGLIBS)
#for optimizing:
#COMPILE.C++ = $(C++) $(OPTFLAGS) -DTIMING -DUSE_GTK
#COMPILE.f = $(FC) $(OPTFLAGS)
#CLASS_LIBS = $(OPTLIBS)

%.o: %.C
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<
all : GUIFourier GUIFiniteDifference GUIExtrapolate GUIQuadrature \
  GUIAdaptiveQuadrature
GUIFourier : GUIFourier.o dffti.o dfftf.o dfftb.o rffti1.o rfftf1.o \
	rfftb1.o radfg.o radf2.o radf3.o radf4.o radf5.o \
	radbg.o radb2.o radb3.o radb4.o radb5.o
	$(C++) -o $@ GUIFourier.o dffti.o dfftf.o dfftb.o \
	radfg.o radf2.o radf3.o radf4.o radf5.o \
	radbg.o radb2.o radb3.o radb4.o radb5.o \
	rffti1.o rfftf1.o rfftb1.o $(CLASS_LIBS) $(GLIBS) \
	-lrfftw -lfftw
	$(MAKE) depends.gnu
GUIFiniteDifference : GUIFiniteDifference.o
	$(C++) -o $@ GUIFiniteDifference.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
GUIExtrapolate : GUIExtrapolate.o
	$(C++) -o $@ GUIExtrapolate.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
GUIQuadrature : GUIQuadrature.o quadrature.o
	$(C++) -o $@ GUIQuadrature.o quadrature.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
GUIAdaptiveQuadrature : GUIAdaptiveQuadrature.o quadrature.o
	$(C++) -o $@ GUIAdaptiveQuadrature.o quadrature.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu

clean : FORCE
	-rm *.o trash* core*
	-rm GUIFourier GUIFiniteDifference GUIExtrapolate GUIQuadrature
	-rm GUIAdaptiveQuadrature
FORCE :

depends.gnu :
	-$(RM) depends.gnu
	cat *.d > depends.gnu
include depends.gnu
