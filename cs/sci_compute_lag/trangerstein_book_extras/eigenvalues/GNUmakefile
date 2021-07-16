FC = gfortran
#FFLAGS = -g -finit-real=INF -fsecond-underscore -Wunused -Werror
FFLAGS = -g -O3
LIBS = -lgfortran
COMPILE.f = $(FC) $(FFLAGS)

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
DEBUGFLAGS = -g -fno-inline -DMEM_DEBUG -DTIMING
OPTFLAGS = -O3 -g -DTIMING
COMPILE.C++ = $(C++) $(DEBUGFLAGS) 
#COMPILE.C++ = $(C++) $(OPTFLAGS)
DGESVD_FILES = dgesvd.o lsame.o xerbla.o ilaenv.o dgelqf.o dorglq.o \
  dgebrd.o dorgbr.o dlange.o dgeqrf.o dlacpy.o dlaset.o dlamch.o dorgqr.o \
  dbdsqr.o dlascl.o dgemm.o dormbr.o dgelq2.o dlarft.o iparmq.o ieeeck.o \
  dlarfb.o dorgl2.o dlabrd.o dgebd2.o dlassq.o dgeqr2.o dorg2r.o dlartg.o \
  dlasr.o dlasv2.o dlas2.o dlasq1.o dormqr.o dormlq.o dlarfg.o dlarf.o \
  drot.o dscal.o dswap.o dgemv.o dtrmv.o disnan.o iladlr.o iladlc.o \
  dcopy.o dtrmm.o dnrm2.o dger.o dlasq2.o dlasrt.o dorm2r.o dorml2.o \
  dlapy2.o dlasq3.o dlaisnan.o dlasq4.o dlasq5.o dlasq6.o daxpy.o ddot.o

CC = gcc -MMD
COMPILE.C = $(CC) -g
#COMPILE.C = $(CC) -O3

%.o: %.C
	$(COMPILE.C++) -c -o $*.o $<
%.o: %.c
	$(COMPILE.C) -c -o $*.o $<
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<
%.o: %.f90
	$(COMPILE.f) -c -o $*.o $<

all: test_dstedc test_dgeev run_dchkhs newton_pinv hw19
test_dstedc : test_dstedc.o 
	$(FC) -o $@ test_dstedc.o -llapack -lblas
test_dgeev : test_dgeev.o dgeev.o dgehrd.o dgehd2.o dhseqr.o dlahqr.o 
	$(FC) -o $@ test_dgeev.o dgeev.o dgehrd.o dgehd2.o dhseqr.o dlahqr.o -llapack -lblas
run_dchkhs : run_dchkhs.o dchkhs.o dlatms.o dlatmr.o dhst01.o dget10.o dget22.o dlatme.o dlafts.o dlasum.o dlatm1.o dlarnd.o dlarot.o dlagge.o dlagsy.o dlatm3.o dlatm2.o dort01.o dlaran.o dlarge.o dlahd2.o
	$(FC) -o $@ run_dchkhs.o dchkhs.o dlatms.o dlatmr.o dhst01.o dget10.o dget22.o dlatme.o dlafts.o dlasum.o dlatm1.o dlarnd.o dlarot.o dlagge.o dlagsy.o dlatm3.o dlatm2.o dort01.o dlaran.o dlarge.o dlahd2.o -llapack -lblas
newton_pinv : newton_pinv.o 
	$(FC) -o $@ newton_pinv.o -lblas
hw19 : hw19.o 
	$(FC) -o $@ hw19.o /usr/lib64/atlas/liblapack.so.3.0
hw_jac : hw_jac.o 
	$(FC) -o $@ hw_jac.o /usr/lib64/atlas/liblapack.so.3.0

clean : FORCE
	-rm *.o core*
	-rm test_dstedc
	-rm test_dstedc test_dgeev run_dchkhs newton_pinv hw19
FORCE :

depends.gnu :
	cat *.d > depends.gnu

include depends.gnu
