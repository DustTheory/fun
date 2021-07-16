#DEBUGFLAGS = -g -fno-implicit-templates -fno-default-inline
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -Wall
FC = gfortran

DEBUGCOMPILEf = $(FC) $(DEBUGFFLAGS)
OPTCOMPILEf = $(FC) $(OPTFFLAGS)

COMPLEX_SOURCES = caxpy.f ccopy.f cdotc.f cdotu.f cgbmv.f cgemm.f cgemv.f \
  cgerc.f cgeru.f chbmv.f chemm.f chemv.f cher.f cher2.f cher2k.f cherk.f \
  chpmv.f chpr.f chpr2.f crotg.f cscal.f csrot.f csscal.f cswap.f csymm.f \
  csyr2k.f csyrk.f ctbmv.f ctbsv.f ctpmv.f ctpsv.f ctrmm.f ctrmv.f \
  ctrsm.f ctrsv.f icamax.f
DOUBLE_REAL_SOURCES = dasum.f daxpy.f dcabs1.f dcopy.f ddot.f dgbmv.f \
  dgemm.f dgemv.f dger.f dnrm2.f drot.f drotg.f drotm.f drotmg.f dsbmv.f \
  dscal.f dsdot.f dspmv.f dspr.f dspr2.f dswap.f dsymm.f dsymv.f dsyr.f \
  dsyr2.f dsyr2k.f dsyrk.f dtbmv.f dtbsv.f dtpmv.f dtpsv.f dtrmm.f dtrmv.f \
  dtrsm.f dtrsv.f dzasum.f dznrm2.f idamax.f
REAL_SOURCES = sasum.f saxpy.f scabs1.f scasum.f scnrm2.f scopy.f sdot.f \
  sdsdot.f sgbmv.f sgemm.f sgemv.f sger.f snrm2.f srot.f srotg.f srotm.f \
  srotmg.f ssbmv.f sscal.f sspmv.f sspr.f sspr2.f sswap.f ssymm.f ssymv.f \
  ssyr.f ssyr2.f ssyr2k.f ssyrk.f stbmv.f stbsv.f stpmv.f stpsv.f strmm.f \
  strmv.f strsm.f strsv.f isamax.f
DOUBLE_COMPLEX_SOURCES = zaxpy.f zcopy.f zdotc.f zdotu.f zdrot.f zdscal.f \
  zgbmv.f zgemm.f zgemv.f zgerc.f zgeru.f zhbmv.f zhemm.f zhemv.f zher.f \
  zher2.f zher2k.f zherk.f zhpmv.f zhpr.f zhpr2.f zrotg.f zscal.f zswap.f \
  zsymm.f zsyr2k.f zsyrk.f ztbmv.f ztbsv.f ztpmv.f ztpsv.f ztrmm.f \
  ztrmv.f ztrsm.f ztrsv.f izamax.f
fSOURCES = lsame.f $(COMPLEX_SOURCES) $(DOUBLE_REAL_SOURCES) \
  $(REAL_SOURCES) $(DOUBLE_COMPLEX_SOURCES)
 
fDEBUGOBJECTS = $(foreach SOURCE, $(fSOURCES), d/$(subst .f,.o,$(SOURCE)))
fOPTOBJECTS = $(foreach SOURCE, $(fSOURCES), o/$(subst .f,.o,$(SOURCE)))

d/%.o : %.f
	$(DEBUGCOMPILEf) -c -o d/$*.o $<
o/%.o : %.f
	$(OPTCOMPILEf) -c -o o/$*.o $<

all :
	-mkdir d
	-mkdir o
	make d/libblas.a
	make o/libblas.a
d/libblas.a : $(fDEBUGOBJECTS) 
	ar -rv d/libblas.a $(fDEBUGOBJECTS)
o/libblas.a : $(fOPTOBJECTS) 
	ar -rv o/libblas.a $(fOPTOBJECTS)

clean :
	-rm d/*.o d/*.a
	-rm o/*.o o/*.a
