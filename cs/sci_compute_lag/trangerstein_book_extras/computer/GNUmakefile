SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))
MINPACK_DIR = $(DIR_BASE)/minpack
ATLAS_DIR = $(DIR_BASE)/ATLAS/MyObj/lib
ATLAS_LIBS = $(ATLAS_DIR)/libf77blas.a $(ATLAS_DIR)/libatlas.a
ATLAS_PTLIBS = $(ATLAS_DIR)/libptf77blas.a $(ATLAS_DIR)/libatlas.a -lpthread

FC = gfortran
FFLAGS = -c -g -finit-real=INF -fsecond-underscore -Wunused -Werror
LIBS = -lgfortran
COMPILE.f = $(FC) $(FFLAGS)

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
DEBUGFLAGS = -g -fno-inline -DMEM_DEBUG -DTIMING
OPTFLAGS = -O3 -g -DTIMING
COMPILE.C++ = $(C++) $(DEBUGFLAGS) 
#COMPILE.C++ = $(C++) $(OPTFLAGS)

CC = gcc -MMD
COMPILE.C = $(CC) -g
#COMPILE.C = $(CC) -O3

%.o: %.C
	$(COMPILE.C++) -c -DDEBUG -o $@ $<
%.o: %.c
	$(COMPILE.C) -c -o $@ $<
%.o: %.f
	$(COMPILE.f) -c -o $@ $<
%.o: %.f90
	$(COMPILE.f) -c -O3 -o $@ $<
all: depends.gnu \
  intOverflowCpp intOverflowC integerOverflowF integerOverflowF90 \
  floatOverflowCpp floatOverflowC realOverflowF realOverflowF90 \
  roundingErrorsCpp roundingErrorsC roundingErrorsF roundingErrorsF90 \
  derivativeCpp derivativeC derivativef derivativef90 \
  quadraticCpp quadraticC quadraticf quadraticf90 \
  exceptionsCpp exceptionsC exceptionsf exceptionsf90 \
  exceptionsCpp_trap exceptionsC_trap exceptionsf_trap exceptionsf90_trap \
  execf execf90 execC execc \
  callingLapackFromFortran callingLapackFromFortran90 callingLapackFromC \
    callingLapackFromCpp \
  callingMinpackFromFortran callingMinpackFromFortran90 \
    callingMinpackFromC callingMinpackFromCpp \
  integrateIO stride copy testMemoryDebugger integrate \
  cio fio cppio \
  bad_integrate fcn fcnOpt multDiv \
  matrixMatrixMultiply_ijk matrixMatrixMultiply_ikj \
    matrixMatrixMultiply_jik matrixMatrixMultiply_jki \
    matrixMatrixMultiply_kij matrixMatrixMultiply_kji \
  matrixMatrixMultiplyDdot_ij matrixMatrixMultiplyDdot_ji \
    matrixMatrixMultiplyDaxpy_jk matrixMatrixMultiplyDaxpy_kj \
  matrixMatrixMultiplyDgemv matrixMatrixMultiplyMatmul \
  matrixMatrixMultiplyDgemmLapack \
  matrixMatrixMultiplyDgemmAtlas matrixMatrixMultiplyDgemmMyAtlas \
    matrixMatrixMultiplyDgemmMyPtAtlas \
  matrixMatrixMultiplyDgemmMyPlasma \
  pthread_hello pthread_dotprod pthread_condvar pthread_barrier \
  depends.gnu

tarfile: FORCE
	tar -cvf tarfile GNUmakefile *.[CHchif] *.f90 *.html *.js

intOverflowCpp : intOverflowCpp.o
	$(C++) -o $@ intOverflowCpp.o 
intOverflowC : intOverflowC.o
	$(CC) -o $@ intOverflowC.o 
integerOverflowF : integerOverflowF.o
	$(FC) -o $@ integerOverflowF.o 
integerOverflowF90 : integerOverflowF90.o
	$(FC) -o $@ integerOverflowF90.o 

floatOverflowCpp : floatOverflowCpp.o
	$(C++) -o $@ floatOverflowCpp.o 
floatOverflowC : floatOverflowC.o
	$(CC) -o $@ floatOverflowC.o 
realOverflowF : realOverflowF.o
	$(FC) -o $@ realOverflowF.o 
realOverflowF90 : realOverflowF90.o
	$(FC) -o $@ realOverflowF90.o 

roundingErrorsCpp : roundingErrorsCpp.o
	$(C++) -o $@ roundingErrorsCpp.o 
roundingErrorsC : roundingErrorsC.o
	$(CC) -o $@ roundingErrorsC.o -lm 
roundingErrorsF : roundingErrorsF.o
	$(FC) -o $@ roundingErrorsF.o
roundingErrorsF90 : roundingErrorsF90.o
	$(FC) -o $@ roundingErrorsF90.o

derivativeCpp : derivativeCpp.o
	$(C++) -o $@ derivativeCpp.o
derivativeC : derivativeC.o
	$(CC) -o $@ derivativeC.o -lm
derivativef : derivativef.o
	$(FC) -o $@ derivativef.o 
derivativef90 : derivativef90.o
	$(FC) -o $@ derivativef90.o 

quadraticCpp : quadraticCpp.o
	$(C++) -o $@ quadraticCpp.o 
quadraticC : quadraticC.o
	$(CC) -o $@ quadraticC.o -lm
quadraticf : quadraticf.o
	$(FC) -o $@ quadraticf.o 
quadraticf90 : quadraticf90.o
	$(FC) -o $@ quadraticf90.o

exceptionsCpp : exceptionsCpp.C
	$(C++) -g -o $@ -D_TRAP_FPE_=0 exceptionsCpp.C 
exceptionsCpp_trap : exceptionsCpp.C
	$(C++) -g -o $@ -D_TRAP_FPE_=1 exceptionsCpp.C 
exceptionsC : exceptionsC.c
	$(CC) -g -o $@ -D_TRAP_FPE_=0 exceptionsC.c -lm 
exceptionsC_trap : exceptionsC.c
	$(CC) -g -o $@ -D_TRAP_FPE_=1 exceptionsC.c -lm 
exceptionsf : exceptionsf.f
	$(FC) -o $@ -g -fno-range-check exceptionsf.f
exceptionsf_trap : exceptionsf.f
	$(FC) -g -fno-range-check -ffpe-trap=invalid,zero,overflow \
	  -o $@ exceptionsf.f
exceptionsf90 : exceptionsf90.f90
	$(FC) -o $@ -g -fno-range-check exceptionsf90.f90
exceptionsf90_trap : exceptionsf90.f90
	$(FC) -g -fno-range-check -ffpe-trap=invalid,zero,overflow \
	  -o $@ exceptionsf90.f90

execf : fmain.o fsub.o
	$(FC) -o $@ fmain.o fsub.o
execf90 : f90sub.o f90main.o 
	$(FC) -o $@ f90main.o f90sub.o 
execC : Cppmain.o Cppsub.o
	$(C++) -o $@ Cppmain.o Cppsub.o
execc : cmain.o csub.o
	$(CC) -o $@ cmain.o csub.o

callingLapackFromFortran : callingLapackFromFortran.o
	$(FC) -o $@ callingLapackFromFortran.o -lblas
callingLapackFromFortran90 : callingLapackFromFortran90.o
	$(FC) -o $@ callingLapackFromFortran90.o -lblas
callingLapackFromC : callingLapackFromC.o
	$(CC) -o $@ callingLapackFromC.o -lblas
callingLapackFromCpp : callingLapackFromCpp.o
	$(C++) -o $@ callingLapackFromCpp.o -lblas

callingMinpackFromFortran : callingMinpackFromFortran.o
	$(FC) -o $@ callingMinpackFromFortran.o $(MINPACK_DIR)/libminpack.a
callingMinpackFromFortran90 : callingMinpackFromFortran90.o
	$(FC) -o $@ callingMinpackFromFortran90.o \
	  $(MINPACK_DIR)/libminpack.a
callingMinpackFromC : callingMinpackFromC.o
	$(CC) -o $@ callingMinpackFromC.o $(MINPACK_DIR)/libminpack.a
callingMinpackFromCpp : callingMinpackFromCpp.o
	$(C++) -o $@ callingMinpackFromCpp.o $(MINPACK_DIR)/libminpack.a

integrateIO : integrate.o integrateIOMain.o TimedObject.o Tracer.o \
  MemoryDebugger.o
	$(C++) -o $@ integrateIOMain.o integrate.o TimedObject.o Tracer.o \
	  MemoryDebugger.o

stride : strideMain.o stride.o TimedObject.o PermutationGenerator.o \
  NumPtr.o MemoryDebugger.o
	$(C++) -o $@ strideMain.o stride.o TimedObject.o \
	PermutationGenerator.o NumPtr.o MemoryDebugger.o

copy : copyMain.o f90copy.o TimedObject.o
	$(C++) -o $@ copyMain.o f90copy.o TimedObject.o -llapack -lblas

testMemoryDebugger : testMemoryDebugger.o MemoryDebugger.o Tracer.o \
testfortroutine.o
	$(C++) -o $@ testMemoryDebugger.o MemoryDebugger.o Tracer.o \
	testfortroutine.o

integrate : mixedMain.o integrate.o Tracer.o TimedObject.o MemoryDebugger.o
	$(C++) -o $@ mixedMain.o integrate.o Tracer.o TimedObject.o \
	  MemoryDebugger.o

cio : cio.o
	$(CC) -o $@ cio.o
fio : fio.o
	$(FC) -o $@ fio.o
cppio : cppio.o
	$(C++) -o $@ cppio.o

bad_integrate : mixedMainBad.o integrate.o Tracer.o TimedObject.o \
  MemoryDebugger.o
	$(C++) -o $@ mixedMainBad.o integrate.o Tracer.o TimedObject.o \
	  MemoryDebugger.o

fcn : fcn.o fcnMain.o TimedObject.o
	$(C++) -o $@ fcnMain.o fcn.o TimedObject.o

fcnOpt : FORCE
	$(FC) -c -g -O3 -o fcn.o fcn.f
	$(C++) -c -g -O3 -o fcnMain.o fcnMain.C
	$(C++) -c -g -O3 -o TimedObject.o TimedObject.C
	$(C++) -o $@ fcnMain.o fcn.o TimedObject.o

multDiv : multDiv.o TimedObject.o
	$(C++) -o $@ multDiv.o TimedObject.o

matrixMatrixMultiply_ijk : matrixMatrixMultiply_ijk.o
	$(FC) -o $@ matrixMatrixMultiply_ijk.o
matrixMatrixMultiply_ikj : matrixMatrixMultiply_ikj.o
	$(FC) -o $@ matrixMatrixMultiply_ikj.o
matrixMatrixMultiply_jik : matrixMatrixMultiply_jik.o
	$(FC) -o $@ matrixMatrixMultiply_jik.o
matrixMatrixMultiply_jki : matrixMatrixMultiply_jki.o
	$(FC) -o $@ matrixMatrixMultiply_jki.o
matrixMatrixMultiply_kij : matrixMatrixMultiply_kij.o
	$(FC) -o $@ matrixMatrixMultiply_kij.o
matrixMatrixMultiply_kji : matrixMatrixMultiply_kji.o
	$(FC) -o $@ matrixMatrixMultiply_kji.o
matrixMatrixMultiplyDdot_ij : matrixMatrixMultiplyDdot_ij.o
	$(FC) -o $@ matrixMatrixMultiplyDdot_ij.o -lblas
matrixMatrixMultiplyDdot_ji : matrixMatrixMultiplyDdot_ji.o
	$(FC) -o $@ matrixMatrixMultiplyDdot_ji.o -lblas
matrixMatrixMultiplyDaxpy_jk : matrixMatrixMultiplyDaxpy_jk.o
	$(FC) -o $@ matrixMatrixMultiplyDaxpy_jk.o -lblas
matrixMatrixMultiplyDaxpy_kj : matrixMatrixMultiplyDaxpy_kj.o
	$(FC) -o $@ matrixMatrixMultiplyDaxpy_kj.o -lblas
matrixMatrixMultiplyDgemv : matrixMatrixMultiplyDgemv.o
	$(FC) -o $@ matrixMatrixMultiplyDgemv.o -lblas
matrixMatrixMultiplyMatmul : matrixMatrixMultiplyMatmul.o
	$(FC) -o $@ matrixMatrixMultiplyMatmul.o
matrixMatrixMultiplyDgemmLapack : matrixMatrixMultiplyDgemm.o
	$(FC) -o $@ matrixMatrixMultiplyDgemm.o -lblas
matrixMatrixMultiplyDgemmAtlas : matrixMatrixMultiplyDgemm.o
	$(FC) -o $@ matrixMatrixMultiplyDgemm.o \
	/usr/lib64/atlas/libf77blas.so.3.0
matrixMatrixMultiplyDgemmMyAtlas : matrixMatrixMultiplyDgemm.o
	$(FC) -o $@ matrixMatrixMultiplyDgemm.o $(ATLAS_LIBS)
matrixMatrixMultiplyDgemmMyPtAtlas : matrixMatrixMultiplyDgemm.o
	$(FC) -o $@ matrixMatrixMultiplyDgemm.o $(ATLAS_PTLIBS)
matrixMatrixMultiplyDgemmMyPlasma : matrixMatrixMultiplyDgemm.o
	$(FC) -o $@ matrixMatrixMultiplyDgemm.o \
	  /usr/lib64/libplasma.so.3.0.0 $(ATLAS_LIBS)

pthread_hello : hello_arg2.o
	$(CC) -o $@ hello_arg2.o -lpthread
pthread_dotprod : dotprod_mutex.o
	$(CC) -o $@ dotprod_mutex.o -lpthread
pthread_condvar : condvar.o
	$(CC) -o $@ condvar.o -lpthread
pthread_barrier : pthread_barrier.o
	$(CC) -o $@  pthread_barrier.o -lpthread

clean : FORCE
	-rm *.o core*
	-rm \
	  intOverflowCpp intOverflowC integerOverflowF integerOverflowF90 \
          floatOverflowCpp floatOverflowC realOverflowF realOverflowF90 \
          roundingErrorsCpp roundingErrorsC roundingErrorsF \
	    roundingErrorsF90 \
	  derivativeCpp derivativeC derivativef derivativef90 \
	  quadraticCpp quadraticC quadraticf quadraticf90 \
	  exceptionsCpp exceptionsC exceptionsf exceptionsf90 \
	  exceptionsCpp_trap exceptionsC_trap exceptionsf_trap \
	    exceptionsf90_trap \
	  execf execf90 execC execc \
	  callingLapackFromFortran callingLapackFromFortran90 \
	    callingLapackFromC callingLapackFromCpp \
	  callingMinpackFromFortran callingMinpackFromFortran90 \
	    callingMinpackFromC callingMinpackFromCpp \
	  integrateIO stride copy testMemoryDebugger integrate \
	  cio fio cppio \
	  bad_integrate fcn fcnOpt multDiv \
	  matrixMatrixMultiply_ijk matrixMatrixMultiply_ikj \
	    matrixMatrixMultiply_jik matrixMatrixMultiply_jki \
	    matrixMatrixMultiply_kij matrixMatrixMultiply_kji \
	  matrixMatrixMultiplyDdot_ij matrixMatrixMultiplyDdot_ji \
	    matrixMatrixMultiplyDaxpy_jk matrixMatrixMultiplyDaxpy_kj \
	  matrixMatrixMultiplyDgemv matrixMatrixMultiplyMatmul \
	  matrixMatrixMultiplyDgemmLapack \
	  matrixMatrixMultiplyDgemmAtlas matrixMatrixMultiplyDgemmMyAtlas \
	    matrixMatrixMultiplyDgemmMyPtAtlas \
	  matrixMatrixMultiplyDgemmMyPlasma
FORCE :

depends.gnu :
	cat *.d > depends.gnu

include depends.gnu
