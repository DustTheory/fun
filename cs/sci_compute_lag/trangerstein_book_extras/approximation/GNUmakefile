SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
LAPACKPP_BASE = $(DIR_BASE)/lapack++
DEBUGLIBS = $(DIR_BASE)/gui/d/libgui.a \
  $(DIR_BASE)/graphics/1d/libgraphics.a \
  $(DIR_BASE)/memdebug/d/libmemdebug.a \
  -lgthread-2.0
LAPACKPP_DEBUG_LIB = $(LAPACKPP_BASE)/d/lapack++.a
OPTLIBS = $(DIR_BASE)/gui/o/libgui.a \
  $(DIR_BASE)/graphics/1o/libgraphics.a \
  $(DIR_BASE)/memdebug/o/libmemdebug.a
LAPACKPP_OPT_LIB = $(LAPACKPP_BASE)/o/lapack++.a
INCLUDE_OPTION = -I$(DIR_BASE)/gui -I$(DIR_BASE)/graphics \
  -I$(DIR_BASE)/memdebug `pkg-config --cflags gtkglext-x11-1.0` \
  -I$(LAPACKPP_BASE)
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
LAPACKPP_LIB = $(LAPACKPP_DEBUG_LIB)

#for optimizing:
#COMPILE.C++ = $(C++) $(OPTFLAGS) -DTIMING -DUSE_GTK
#COMPILE.f = $(FC) $(OPTFLAGS)
#CLASS_LIBS = $(OPTLIBS)

%.o: %.C
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<

all : erf interpolate GUIInterpolate piecewisePolynomial \
  GUIPiecewisePolynomial bSpline tensionSpline \
  fourier GUIFourier testBinomialCoefficient \
  testDaubechiesFilterBank testSplineFilterBank testSplineScalingFunction \
  testDaubechiesScalingFunction testDaubechiesFilterBank2 \
  testSplineFilterBank2 GUIWavelet
erf : erf.o
	$(COMPILE.C++) -o $@ erf.o $(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
interpolate : interpolate.o divided_difference.o
	$(COMPILE.C++) -o $@ interpolate.o divided_difference.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
GUIInterpolate: GUIInterpolate.o divided_difference.o
	$(COMPILE.C++) -o $@ GUIInterpolate.o divided_difference.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
piecewisePolynomial : piecewise_poly.o piecewisePolynomial.o
	$(COMPILE.C++) -o $@ piecewisePolynomial.o piecewise_poly.o \
        $(DIR_BASE)/graphics/1d/libgraphics.a \
        $(DIR_BASE)/memdebug/d/libmemdebug.a \
        -lgthread-2.0 \
	$(GLIBS) $(FLIBS) -llapack
	$(MAKE) depends.gnu
GUIPiecewisePolynomial : GUIPiecewisePolynomial.o piecewise_poly.o
	$(COMPILE.C++) -o $@ GUIPiecewisePolynomial.o piecewise_poly.o \
	$(CLASS_LIBS) $(GLIBS) $(FLIBS) -llapack
	$(MAKE) depends.gnu
bSpline : bSpline.o
	$(COMPILE.C++) -o $@ bSpline.o $(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
tensionSpline : tensionSpline.o
	$(COMPILE.C++) -o $@ tensionSpline.o \
	$(CLASS_LIBS) $(GLIBS) -llapack
	$(MAKE) depends.gnu
#hermiteCubicSpline : hermiteCubicSpline.o lr.o spline.o
#	$(COMPILE.C++) -o $@ hermiteCubicSpline.o lr.o spline.o \
#	$(CLASS_LIBS) $(GLIBS)
#	$(MAKE) depends.gnu
#smoothCubicSpline : smoothCubicSpline.o
#	$(COMPILE.C++) -o $@ smoothCubicSpline.o \
#	$(CLASS_LIBS) $(GLIBS)
#	$(MAKE) depends.gnu
fourier : fourier.o dffti.o dfftf.o dfftb.o rffti1.o rfftf1.o rfftb1.o \
	radfg.o radf2.o radf3.o radf4.o radf5.o \
	radbg.o radb2.o radb3.o radb4.o radb5.o
	$(COMPILE.C++) -o $@ fourier.o dffti.o dfftf.o dfftb.o \
	radfg.o radf2.o radf3.o radf4.o radf5.o \
	radbg.o radb2.o radb3.o radb4.o radb5.o \
	rffti1.o rfftf1.o rfftb1.o $(CLASS_LIBS) $(GLIBS) \
	-lrfftw -lfftw
	$(MAKE) depends.gnu
GUIFourier : GUIFourier.o dffti.o dfftf.o dfftb.o rffti1.o rfftf1.o \
	rfftb1.o radfg.o radf2.o radf3.o radf4.o radf5.o \
	radbg.o radb2.o radb3.o radb4.o radb5.o
	$(COMPILE.C++) -o $@ GUIFourier.o dffti.o dfftf.o dfftb.o \
	radfg.o radf2.o radf3.o radf4.o radf5.o \
	radbg.o radb2.o radb3.o radb4.o radb5.o \
	rffti1.o rfftf1.o rfftb1.o $(CLASS_LIBS) $(GLIBS) \
	-lrfftw -lfftw
	$(MAKE) depends.gnu
testBinomialCoefficient : testBinomialCoefficient.o BinomialCoefficient.o
	$(COMPILE.C++) -o $@ testBinomialCoefficient.o \
	BinomialCoefficient.o $(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
testDaubechiesFilterBank : BinomialCoefficient.o Signal.o Filter.o \
	  DaubechiesFilterBank.o testDaubechiesFilterBank.o
	$(COMPILE.C++) -o $@ testDaubechiesFilterBank.o \
	  BinomialCoefficient.o Signal.o Filter.o DaubechiesFilterBank.o \
	$(CLASS_LIBS) $(GLIBS) $(LAPACKPP_LIB) -llapack -lblas
	$(MAKE) depends.gnu
testSplineFilterBank : BinomialCoefficient.o Signal.o Filter.o \
	  SplineFilterBank.o testSplineFilterBank.o
	$(COMPILE.C++) -o $@ testSplineFilterBank.o BinomialCoefficient.o \
	  Signal.o Filter.o SplineFilterBank.o \
	$(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
testSplineScalingFunction : testSplineScalingFunction.o \
	  BinomialCoefficient.o Signal.o Filter.o ScalingFunction.o \
	  SplineFilterBank.o Wavelet.o 
	$(COMPILE.C++) -o $@ testSplineScalingFunction.o \
	  BinomialCoefficient.o Signal.o Filter.o ScalingFunction.o \
	  SplineFilterBank.o Wavelet.o \
	$(CLASS_LIBS) $(GLIBS) $(LAPACKPP_LIB) -llapack -lblas
	$(MAKE) depends.gnu
testDaubechiesScalingFunction : testDaubechiesScalingFunction.o \
	  BinomialCoefficient.o Signal.o Filter.o ScalingFunction.o \
	  DaubechiesFilterBank.o Wavelet.o
	$(COMPILE.C++) -o $@ testDaubechiesScalingFunction.o \
	  BinomialCoefficient.o Signal.o Filter.o ScalingFunction.o \
	  DaubechiesFilterBank.o Wavelet.o \
	$(CLASS_LIBS) $(GLIBS) $(LAPACKPP_LIB) -llapack -lblas
	$(MAKE) depends.gnu
testDaubechiesFilterBank2 : BinomialCoefficient.o Signal.o Filter.o \
	  DaubechiesFilterBank.o ScalingFunction.o \
	  testDaubechiesFilterBank2.o
	$(COMPILE.C++) -o $@ testDaubechiesFilterBank2.o \
	  BinomialCoefficient.o Signal.o Filter.o DaubechiesFilterBank.o \
	  ScalingFunction.o $(CLASS_LIBS) $(GLIBS) $(LAPACKPP_LIB) \
	  -llapack -lblas
	$(MAKE) depends.gnu
testSplineFilterBank2 : BinomialCoefficient.o Signal.o Filter.o \
	  SplineFilterBank.o ScalingFunction.o testSplineFilterBank2.o
	$(COMPILE.C++) -o $@ testSplineFilterBank2.o \
	  BinomialCoefficient.o Signal.o Filter.o SplineFilterBank.o \
	  ScalingFunction.o $(CLASS_LIBS) $(GLIBS) $(LAPACKPP_LIB) \
	  -llapack -lblas
	$(MAKE) depends.gnu
GUIWavelet : GUIWavelet.o
	$(COMPILE.C++) -o $@ GUIWavelet.o $(CLASS_LIBS) $(GLIBS)
	$(MAKE) depends.gnu
clean : FORCE
	-rm *.o trash* core*
	-rm erf interpolate GUIInterpolate piecewisePolynomial
	-rm GUIPiecewisePolynomial bSpline tensionSpline 
	-rm fourier GUIFourier testBinomialCoefficient
	-rm testDaubechiesFilterBank testSplineFilterBank 
	-rm testSplineScalingFunction testDaubechiesScalingFunction 
	-rm testDaubechiesFilterBank2 testSplineFilterBank2 GUIWavelet
FORCE :

depends.gnu :
	-$(RM) depends.gnu
	cat *.d > depends.gnu
include depends.gnu
