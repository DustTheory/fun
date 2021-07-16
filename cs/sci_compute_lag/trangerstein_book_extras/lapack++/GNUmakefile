#"%W%\t%G%"
#SD is the directory where the source code is to be found
SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
CLASSDIR = $(DIR_BASE)
DEBUGLIBS = $(CLASSDIR)/memdebug/d/libmemdebug.a
OPTLIBS = $(CLASSDIR)/memdebug/o/libmemdebug.a
INCLUDE_OPTION = -I$(DIR_BASE)/memdebug

#DEBUGFLAGS = -g -fno-implicit-templates -fno-default-inline
DEBUGFLAGS = -g -DDEBUG -DMEM_DEBUG
OPTFLAGS = -O3
CPLUSPLUS = g++ -Wnon-virtual-dtor -Woverloaded-virtual -g -MMD
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -fsecond-underscore -Wall
CC = gcc -g -MMD
FC = gfortran
LINKC = $(CPLUSPLUS)

DEBUGCOMPILEC = $(CPLUSPLUS) $(DEBUGFLAGS)
DEBUGCOMPILEc = $(CC) $(DEBUGFLAGS)
DEBUGCOMPILEf = $(FC) $(DEBUGFFLAGS)
OPTCOMPILEC = $(CPLUSPLUS) $(OPTFLAGS) 
OPTCOMPILEc = $(CC) $(OPTFLAGS) 
OPTCOMPILEf = $(FC) $(OPTFFLAGS)

#fSOURCES = cray.f
#fSOURCES =
fSOURCES = cgbamv.f cgbconnp.f cgbtf2np.f cgbtrsnp.f cgeamv.f \
  cgtconnp.f cgtmv.f cgtrfsrnp.f cgtrfsr.f cgtrfsnp.f cgtsvnp.f \
  cgtsvr.f cgtsvrnp.f cgttrfnp.f cgttrsnp.f cgtts2np.f chbamv.f chtmv.f \
  cla_heamv.f cla_porpvgrw.f cla_syrpvgrw.f cptrfsr.f dgbamv.f \
  dgbconnp.f dgbtf2np.f dgbtrsnp.f dgtconnp.f dgtmv.f dgtrfsnp.f \
  dgtrfsr.f dgtrfsrnp.f dgtsvnp.f dgtsvr.f dgtsvrnp.f dgttrfnp.f \
  dgttrsnp.f dgtts2np.f dla_geamv.f dla_porpvgrw.f dla_syamv.f \
  dla_syrpvgrw.f dptrfsr.f dsbamv.f dstmv.f idmax.f idmin.f idsumn.f \
  idsump.f ismin.f issumn.f issump.f sgbamv.f sgbconnp.f sgbtf2np.f \
  sgbtrsnp.f sgtconnp.f sgtmv.f sgtrfsnp.f sgtrfsr.f sgtrfsrnp.f \
  sgtsvnp.f sgtsvr.f sgtsvrnp.f sgttrfnp.f sgttrsnp.f sla_geamv.f \
  sla_porpvgrw.f sla_syamv.f sla_syrpvgrw.f sptrfsr.f ssbamv.f sstmv.f \
  zgtmv.f zhtmv.f zgbamv.f zgbconnp.f zgbtf2np.f zgbtrsnp.f zgeamv.f \
  zgtconnp.f zgtrfsnp.f zgtrfsr.f zgtrfsrnp.f zgtsvnp.f zgtsvr.f \
  zgtsvrnp.f zgttrfnp.f zgttrsnp.f zgtts2np.f zhbamv.f zla_heamv.f \
  zla_porpvgrw.f zla_syrpvgrw.f zptrfsr.f
cSOURCES = 
#cSOURCES = cgbamv.c cgbconnp.c cgbtf2np.c cgbtrsnp.c cgeamv.c \
#  cgtconnp.c cgtmv.c cgtrfsrnp.c cgtrfsr.c cgtrfsnp.c cgtsvnp.c \
#  cgtsvr.c cgtsvrnp.c cgttrfnp.c cgttrsnp.c cgtts2np.c chbamv.c chtmv.c \
#  cla_heamv.c cla_porpvgrw.c cla_syrpvgrw.c cptrfsr.c dgbamv.c \
#  dgbconnp.c dgbtf2np.c dgbtrsnp.c dgtconnp.c dgtmv.c dgtrfsnp.c \
#  dgtrfsr.c dgtrfsrnp.c dgtsvnp.c dgtsvr.c dgtsvrnp.c dgttrfnp.c \
#  dgttrsnp.c dgtts2np.c dla_geamv.c dla_porpvgrw.c dla_syamv.c \
#  dla_syrpvgrw.c dptrfsr.c dsbamv.c dstmv.c idmax.c idmin.c idsumn.c \
#  idsump.c ismin.c issumn.c issump.c sgbamv.c sgbconnp.c sgbtf2np.c \
#  sgbtrsnp.c sgtconnp.c sgtmv.c sgtrfsnp.c sgtrfsr.c sgtrfsrnp.c \
#  sgtsvnp.c sgtsvr.c sgtsvrnp.c sgttrfnp.c sgttrsnp.c sla_geamv.c \
#  sla_porpvgrw.c sla_syamv.c sla_syrpvgrw.c sptrfsr.c ssbamv.c sstmv.c \
#  zgtmv.c zhtmv.c zgbamv.c zgbconnp.c zgbtf2np.c zgbtrsnp.c zgeamv.c \
#  zgtconnp.c zgtrfsnp.c zgtrfsr.c zgtrfsrnp.c zgtsvnp.c zgtsvr.c \
#  zgtsvrnp.c zgttrfnp.c zgttrsnp.c zgtts2np.c zhbamv.c zla_heamv.c \
#  zla_porpvgrw.c zla_syrpvgrw.c zptrfsr.c
CSOURCES = Vector.C Matrix.C SquareMatrix.C \
  TrapezoidalMatrix.C \
  OrthogonalMatrix.C SymmetricMatrix.C \
  BandMatrix.C SpecializedMatrix.C \
  GaussianFactorization.C CholeskyFactorization.C \
  GramSchmidtQRFactorization.C HouseholderQRFactorization.C \
  CompleteOrthogonalDecomposition.C SingularValueDecomposition.C \
  LaDouble.C LaDoubleComplex.C LaFloat.C LaComplex.C
# LaEnum.C LinearProgram.C \
 
CDEBUGOBJECTS = $(foreach SOURCE, $(CSOURCES), d/$(subst .C,.o,$(SOURCE)))
ifneq ( ,$(cSOURCES))
  cDEBUGOBJECTS = $(foreach SOURCE, $(cSOURCES), d/$(subst .c,.o,$(SOURCE)))
endif
ifneq ( ,$(fSOURCES))
  fDEBUGOBJECTS = $(foreach SOURCE, $(fSOURCES), d/$(subst .f,.o,$(SOURCE)))
endif
COPTOBJECTS = $(foreach SOURCE, $(CSOURCES), o/$(subst .C,.o,$(SOURCE)))
ifneq ( ,$(cSOURCES))
  cOPTOBJECTS = $(foreach SOURCE, $(cSOURCES), o/$(subst .c,.o,$(SOURCE)))
endif
ifneq ( ,$(fSOURCES))
  fOPTOBJECTS = $(foreach SOURCE, $(fSOURCES), o/$(subst .f,.o,$(SOURCE)))
else
  DEBUGLIBS := $(DEBUGLIBS) -lf2c
  OPTLIBS := $(OPTLIBS) -lf2c
endif

d/%.o: %.C
	$(DEBUGCOMPILEC) $(INCLUDE_OPTION) -c -o d/$*.o $< 
d/%.o: %.f
	$(DEBUGCOMPILEf) -c -o d/$*.o $<
d/%.o: %.c
	$(DEBUGCOMPILEc) -c -o d/$*.o $<
o/%.o: %.C
	$(OPTCOMPILEC) $(INCLUDE_OPTION) -c -o o/$*.o $< 
o/%.o: %.f
	$(OPTCOMPILEf) -c -o o/$*.o $<
o/%.o: %.c
	$(OPTCOMPILEc) -c -o o/$*.o $<

all : d/lapack++.a o/lapack++.a
	@echo cSOURCES = $(cSOURCES)
	@echo fSOURCES = $(fSOURCES)
	-rm depends.gnu
	-cat d/*.d o/*.d > depends.gnu
d/lapack++.a : $(CDEBUGOBJECTS) $(cDEBUGOBJECTS) $(fDEBUGOBJECTS) 
	ar -rv d/lapack++.a $(CDEBUGOBJECTS) $(cDEBUGOBJECTS) \
	  $(fDEBUGOBJECTS)
o/lapack++.a : $(COPTOBJECTS) $(cOPTOBJECTS) $(fOPTOBJECTS) 
	ar -rv o/lapack++.a $(COPTOBJECTS) $(cOPTOBJECTS) $(fOPTOBJECTS)

latestd : d/testmain_double.o d/lapack++.a
	  $(LINKC) -o $@ d/testmain_double.o \
	  d/lapack++.a $(DEBUGLIBS) -llapack -lblas
latests : d/testmain_float.o d/lapack++.a
	  $(LINKC) -o $@ d/testmain_float.o d/lapack++.a $(DEBUGLIBS) \
	    -llapack -lblas

simplex: d/simplex.o
	$(LINKC) -o $@ d/simplex.o d/lapack++.a $(DEBUGLIBS)
dual_simplex: d/dual_simplex.o
	$(LINKC) -o $@ d/dual_simplex.o d/lapack++.a $(DEBUGLIBS)

clean :
	-rm d/*.o d/*.a d/*.d
	-rm o/*.o o/*.a o/*.d
zip :
	tar -cvf tarfile *.[CHchf] GNUmakefile
	gzip tarfile
depends.gnu :
	touch depends.gnu
	-mkdir d
	-mkdir o

include depends.gnu
