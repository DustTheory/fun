SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
DEBUGLIBS = $(DIR_BASE)/memdebug/d/libmemdebug.a
OPTLIBS = $(DIR_BASE)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(DIR_BASE)/memdebug -I$(DIR_BASE)/lapack++
#ATLAS_PTLIBS = /usr/lib64/atlas/liblapack.so.3 /usr/lib64/atlas/libptf77blas.so.3
#  /usr/lib64/atlas/libatlas.so.3 -lpthread /usr/lib64/libgfortran.so.3
ATLAS_DIR = $(DIR_BASE)/ATLAS/MyObj/lib
ATLAS_PTLIBS = $(ATLAS_DIR)/libptlapack.a $(ATLAS_DIR)/libf77blas.a $(ATLAS_DIR)/libcblas.a $(ATLAS_DIR)/libatlas.a -lpthread -lgfortran -llapack

FC = gfortran
#FFLAGS = -c -g -fimplicit-none -fno-f2c -finit-real=INF -fno-underscoring -Wunused -Werror
FFLAGS = -c -g -fimplicit-none -fno-f2c -finit-real=INF -fsecond-underscore -Wunused -Werror
LIBS = -lgfortran
COMPILE.f = $(FC) $(FFLAGS)

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
#C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD
DEBUGFLAGS = -g 
OPTFLAGS = -O3 -g 
COMPILE.C++ = $(C++) $(DEBUGFLAGS)

CC = gcc -MMD
COMPILE.C = $(CC) $(DEBUGFLAGS)

GaussianCaseStudy.o: GaussianCaseStudy.C
	$(C++) $(OPTFLAGS) -c -o $*.o $<
%.o: %.c
	$(COMPILE.C) -c -o $*.o $<
%.o: %.C
	$(COMPILE.C++) -c -o $*.o $<
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<
%.o: %.f90
	$(COMPILE.f) -c -o $*.o $<

all : forwardSubF forwardSubF90 forwardSubCpp forwardSubC \
  backSubF backSubF90 backSubCpp backSubC \
  gaussElimFullPivF gaussElimFullPivF90 gaussElimFullPivC \
    gaussElimFullPivCpp \
  gaussElimPartialPivF gaussElimPartialPivF90 gaussElimPartialPivC \
    gaussElimPartialPivCpp \
  GaussianCaseStudy \
  gaussScaleImproveF gaussScaleImproveF90 gaussScaleImproveC \
    gaussScaleImproveCpp \
  ImprovementCaseStudy testVector hw5

# dblat1 dblat2 dblat3 testLaPack0 testLaPack1 testLaPack2 \
# testLaPack3 testLaPack solveNoPiv solveFullPiv solvePartialPiv \
# testVector testMatrix solveMatrix hw4 hw5

forwardSubF : forwardSubF.o
	$(FC) -o $@ forwardSubF.o -lblas
forwardSubF90 : forwardSubF90.o
	$(FC) -o $@ forwardSubF90.o -lblas
forwardSubCpp : forwardSubCpp.o 
	$(C++) -o $@ forwardSubCpp.o -lblas
forwardSubC : forwardSubC.o 
	$(CC) -o $@ forwardSubC.o -lblas

backSubF : backSubF.o
	$(FC) -o $@ backSubF.o -lblas
backSubF90 : backSubF90.o
	$(FC) -o $@ backSubF90.o -lblas
backSubCpp : backSubCpp.o 
	$(C++) -o $@ backSubCpp.o -lblas
backSubC : backSubC.o 
	$(CC) -o $@ backSubC.o -lblas

gaussElimFullPivF : gaussElimFullPivF.o 
	$(FC) -o $@ gaussElimFullPivF.o -llapack -lblas
gaussElimFullPivF90 : gaussElimFullPivF90.o 
	$(FC) -o $@ gaussElimFullPivF90.o -llapack -lblas
gaussElimFullPivC : gaussElimFullPivC.o 
	$(CC) -o $@ gaussElimFullPivC.o -llapack -lblas
gaussElimFullPivCpp : gaussElimFullPivCpp.o 
	$(C++) -o $@ gaussElimFullPivCpp.o $(DEBUGLIBS) \
	  $(DIR_BASE)/lapack++/d/lapack++.a -llapack -lblas

gaussElimPartialPivF : gaussElimPartialPivF.o 
	$(FC) -o $@ gaussElimPartialPivF.o -llapack -lblas
gaussElimPartialPivF90 : gaussElimPartialPivF90.o 
	$(FC) -o $@ gaussElimPartialPivF90.o -llapack -lblas
gaussElimPartialPivC : gaussElimPartialPivC.o 
	$(CC) -o $@ gaussElimPartialPivC.o -llapack -lblas
gaussElimPartialPivCpp : gaussElimPartialPivCpp.o 
	$(C++) -o $@ gaussElimPartialPivCpp.o $(DEBUGLIBS) \
	  $(DIR_BASE)/lapack++/d/lapack++.a -llapack -lblas

GaussianCaseStudy : GaussianCaseStudy.o 
	$(C++) -o $@ GaussianCaseStudy.o $(OPTLIBS) \
	  $(DIR_BASE)/lapack++/o/lapack++.a $(ATLAS_PTLIBS)

gaussScaleImproveF : gaussScaleImproveF.o 
	$(FC) -o $@ gaussScaleImproveF.o -llapack -lblas
gaussScaleImproveF90 : gaussScaleImproveF90.o 
	$(FC) -o $@ gaussScaleImproveF90.o -llapack -lblas
gaussScaleImproveC : gaussScaleImproveC.o 
	$(CC) -o $@ gaussScaleImproveC.o -llapack -lblas
gaussScaleImproveCpp : gaussScaleImproveCpp.o 
	$(C++) -o $@ gaussScaleImproveCpp.o $(DEBUGLIBS) \
	  $(DIR_BASE)/lapack++/d/lapack++.a -llapack -lblas

ImprovementCaseStudy : ImprovementCaseStudy.o 
	$(C++) -o $@ ImprovementCaseStudy.o $(OPTLIBS) \
	  $(DIR_BASE)/lapack++/o/lapack++.a $(ATLAS_PTLIBS)

#dblat1 : dblat1.o ddot.o daxpy.o drotg.o drot.o dcopy.o dswap.o dnrm2.o dasum.o dscal.o idamax.o
#	$(FC) -o $@ dblat1.o ddot.o daxpy.o drotg.o drot.o dcopy.o dswap.o dnrm2.o dasum.o dscal.o idamax.o -lblas
#dblat2 : dblat2.o
#	$(FC) -o $@ dblat2.o -lblas
#dblat3 : dblat3.o
#	$(FC) -o $@ dblat3.o -lblas
#testLaPack0 : testLaPack0.o 
#	$(FC) -o $@ testLaPack0.o -llapack -lblas
#testLaPack1 : testLaPack1.o
#	$(FC) -o $@ testLaPack1.o -llapack -lblas
#testLaPack2 : testLaPack2.o dgelq2.o dgeql2.o dgeqr2.o dgerq2.o
#	$(FC) -o $@ testLaPack2.o dgelq2.o dgeql2.o dgeqr2.o dgerq2.o -llapack -lblas
#testLaPack3 : testLaPack3.o dgelq2.o dgeql2.o dgeqr2.o dgerq2.o
#	$(FC) -o $@ testLaPack3.o dgelq2.o dgeql2.o dgeqr2.o dgerq2.o -llapack -lblas
#testLaPack : testLaPack.o
#	$(FC) -o $@ testLaPack.o -llapack
#solveNoPiv : solveNoPiv.o blas1.o forwardSubWithBlas.o \
#  backSubWithBlas.o gaussElimWithBlas.o
#	$(C++) -o $@ solveNoPiv.o blas1.o forwardSubWithBlas.o \
#	backSubWithBlas.o gaussElimWithBlas.o $(DEBUGLIBS)
#	$(MAKE) depends.gnu
#solveFullPiv : solveFullPiv.o blas1.o forwardSubWithBlas.o \
#  backSubWithBlas.o gaussElimFullPiv.o sswap.o
#	$(C++) -o $@ solveFullPiv.o blas1.o forwardSubWithBlas.o \
#	backSubWithBlas.o gaussElimFullPiv.o sswap.o $(DEBUGLIBS)
#	$(MAKE) depends.gnu
#solvePartialPiv : solvePartialPiv.o blas1.o sgetrf_simple.o \
#  sgetrs_simple.o forwardSubWithBlas.o backSubWithBlas.o sswap.o
#	$(C++) -o $@ solvePartialPiv.o blas1.o sgetrf_simple.o \
#	sgetrs_simple.o forwardSubWithBlas.o backSubWithBlas.o sswap.o \
#	$(DEBUGLIBS)
#	$(MAKE) depends.gnu
testVector : testVector.o vector.o
	$(C++) -o $@ testVector.o vector.o
	$(MAKE) depends.gnu
#testMatrix : testMatrix.o matrix.o vector.o strsv_simple.o \
#  gaussElimFullPiv.o blas1.o sswap.o
#	$(C++) -o $@ testMatrix.o matrix.o vector.o strsv_simple.o \
#	gaussElimFullPiv.o blas1.o sswap.o $(DEBUGLIBS)
#	$(MAKE) depends.gnu
#solveMatrix : solveMatrix.o matrix.o vector.o strsv_simple.o \
#  gaussElimFullPiv.o blas1.o sswap.o
#	$(C++) -o $@ solveMatrix.o matrix.o vector.o strsv_simple.o \
#	gaussElimFullPiv.o blas1.o sswap.o $(DEBUGLIBS)
#	$(MAKE) depends.gnu
hw5 : hw5.o
	$(C++) -o $@ hw5.o $(DIR_BASE)/lapack++/d/lapack++.a\
	  $(DEBUGLIBS) -llapack -lblas
clean : FORCE
	-rm *.o trash* core*
	-rm forwardSubF forwardSubF90 forwardSubCpp forwardSubC
	-rm backSubF backSubF90 backSubCpp backSubC 
	-rm gaussElimFullPivF gaussElimFullPivF90 gaussElimFullPivC \
	  gaussElimFullPivCpp
	-rm gaussElimPartialPivF gaussElimPartialPivF90 \
	  gaussElimPartialPivC gaussElimPartialPivCpp
	-rm GaussianCaseStudy
	-rm gaussScaleImproveF gaussScaleImproveF90 gaussScaleImproveC \
	  gaussScaleImproveCpp \
	-rm pthread_condvar pthread_dotprod pthread_hello \
	-rm ImprovementCaseStudy testVector hw5
FORCE :

depends.gnu :
	-$(RM) depends.gnu
	cat *.d > depends.gnu

include depends.gnu
