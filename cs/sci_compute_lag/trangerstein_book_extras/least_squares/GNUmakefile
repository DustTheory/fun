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
DEBUGFLAGS = -g -DDEBUG -DMEM_DEBUG 
OPTFLAGS = -O3 -g 
COMPILE.C++ = $(C++) $(DEBUGFLAGS)

CC = gcc -MMD
COMPILE.C = $(CC) $(DEBUGFLAGS)

all : LeastSquaresCaseStudy depends.gnu

LeastSquaresCaseStudy.o: LeastSquaresCaseStudy.C
	$(C++) $(OPTFLAGS) -c -o $*.o $<
%.o: %.c
	$(COMPILE.C) -c -o $*.o $<
%.o: %.C
	$(COMPILE.C++) -c -o $*.o $<
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<
%.o: %.f90
	$(COMPILE.f) -c -o $*.o $<

LeastSquaresCaseStudy : LeastSquaresCaseStudy.o
	$(C++) -o $@ LeastSquaresCaseStudy.o $(OPTLIBS) \
	  $(DIR_BASE)/lapack++/o/lapack++.a $(ATLAS_PTLIBS)
clean : FORCE
	-rm *.o trash* core*
	-rm LeastSquaresCaseStudy
FORCE :

depends.gnu :
	-$(RM) depends.gnu
	cat *.d > depends.gnu

include depends.gnu
