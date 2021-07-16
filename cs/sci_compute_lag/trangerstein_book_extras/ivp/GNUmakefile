SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
LAPACKPP_BASE = $(DIR_BASE)/lapack++

DEBUGLIBS = $(DIR_BASE)/gui/d/libgui.a \
  $(DIR_BASE)/graphics/1d/libgraphics.a \
  $(DIR_BASE)/memdebug/d/libmemdebug.a \
  -lgthread-2.0
OPTLIBS = $(DIR_BASE)/gui/o/libgui.a \
  $(DIR_BASE)/graphics/1o/libgraphics.a \
  $(DIR_BASE)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(SD) -I$(DIR_BASE)/gui -I$(DIR_BASE)/graphics \
  -I$(DIR_BASE)/memdebug -I$(LAPACKPP_BASE)\
  `pkg-config --cflags gtkglext-x11-1.0`
GLIBS = `pkg-config --libs gtkglext-x11-1.0`

#DEAL_BASE = /home/johnt/deal_new

FC = gfortran
DEBUGFFLAGS = -c -g -fimplicit-none -fno-f2c -finit-real=INF -fsecond-underscore   -Werror
FLIBS = -lgfortran

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
DEBUGFLAGS = -g
OPTFLAGS = -O3 -g

CC = gcc -MMD

#for debugging:
COMPILE.C++ = $(C++) $(DEBUGFLAGS) -DTIMING -DMEM_DEBUG -DDEBUG -DUSE_GTK -D    NO_THREAD
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
all : GUIAbsoluteStability GUILinearMultistep GUIVode GUICVode \
  GUIDeferredCorrection GUIDCAbsoluteStability GUIRKAbsoluteStability \
  GUINonlinearStability
GUIAbsoluteStability : GUIAbsoluteStability.o BinomialCoefficient.o \
	  integrate.o
	$(C++) -o $@ GUIAbsoluteStability.o BinomialCoefficient.o \
	  integrate.o \
	$(CLASS_LIBS) $(LAPACKPP_LIB) $(GLIBS) $(FLIBS) -llapack -lblas
	$(MAKE) depends.gnu

GUILinearMultistep : GUILinearMultistep.o integrate.o
	$(C++) -o $@ GUILinearMultistep.o integrate.o \
	$(CLASS_LIBS) $(GLIBS) $(FLIBS)
	$(MAKE) depends.gnu

GUIVode : GUIVode.o vode.o dgefa.o dgbfa.o dgesl.o dgbsl.o
	$(C++) -o $@ GUIVode.o vode.o dgefa.o dgbfa.o dgesl.o dgbsl.o \
	$(CLASS_LIBS) $(GLIBS) $(FLIBS) -lblas
	$(MAKE) depends.gnu

GUICVode : GUICVode.o band.o cvband.o cvbandpre.o cvdense.o cvdiag.o \
  cvode.o cvspgmr.o dense.o iterativ.o llnlmath.o nvector.o spgmr.o 
	$(C++) -o $@ GUICVode.o band.o cvband.o cvbandpre.o cvdense.o \
	cvdiag.o cvode.o cvspgmr.o dense.o iterativ.o llnlmath.o \
	nvector.o spgmr.o \
	$(CLASS_LIBS) $(GLIBS) 
	$(MAKE) depends.gnu

GUIDeferredCorrection : GUIDeferredCorrection.o NumPtr.o Polynomial.o \
  Vector.o LaDouble.o Tensor.o Quadrature.o deferred_correction.o
	$(C++) -o $@ GUIDeferredCorrection.o NumPtr.o Polynomial.o \
	Vector.o LaDouble.o Tensor.o Quadrature.o deferred_correction.o \
	$(CLASS_LIBS) $(GLIBS) -lblas -llapack
	$(MAKE) depends.gnu

GUIDCAbsoluteStability : GUIDCAbsoluteStability.o Polynomial.o \
  Quadrature.o deferred_correction.o dogleg.o dpmpar.o enorm.o fdjac1.o \
  hybrd.o qform.o qrfac.o r1mpyq.o r1updt.o
	$(C++) -o $@ GUIDCAbsoluteStability.o Polynomial.o Quadrature.o \
	deferred_correction.o dogleg.o dpmpar.o enorm.o fdjac1.o hybrd.o \
	qform.o qrfac.o r1mpyq.o r1updt.o \
	$(CLASS_LIBS) $(LIBS) $(GLIBS) -llapack
	$(MAKE) depends.gnu

GUIRKAbsoluteStability : GUIRKAbsoluteStability.o deferred_correction.o
	$(C++) -o $@ GUIRKAbsoluteStability.o \
	deferred_correction.o \
	$(CLASS_LIBS) $(LIBS) $(GLIBS)
	$(MAKE) depends.gnu

GUINonlinearStability : GUINonlinearStability.o
	$(C++) -o $@ GUINonlinearStability.o \
	$(CLASS_LIBS) $(LIBS) $(GLIBS) 
	$(MAKE) depends.gnu

clean : FORCE
	-rm *.o trash* core*
	-rm GUIAbsoluteStability GUILinearMultistep GUIVode GUICVode
	-rm GUIDeferredCorrection GUIDCAbsoluteStability 
	-rm GUIRKAbsoluteStability GUINonlinearStability
FORCE :

depends.gnu :
	-$(RM) depends.gnu
	cat *.d > depends.gnu
include depends.gnu
