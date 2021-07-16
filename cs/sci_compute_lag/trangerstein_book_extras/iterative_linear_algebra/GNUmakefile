SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))

FC = gfortran
FFLAGS = -g -fimplicit-none -fno-f2c -finit-real=INF -fsecond-underscore -Wunused -Werror
LIBS = -lgfortran
COMPILE.f = $(FC) $(FFLAGS)

C++ = g++ -Wnon-virtual-dtor -Woverloaded-virtual -MMD $(INCLUDE_OPTION)
DEBUGFLAGS = -g -fno-inline -DDEBUG -DINDEF -DMEM_DEBUG -DUSE_GTK -DNO_THREAD -DTIMING
OPTFLAGS = -O3 -g -DUSE_GTK -DNO_THREAD -DTIMING
COMPILE.C++ = $(C++) $(DEBUGFLAGS)

CC = gcc -MMD
COMPILE.C = $(CC) -g

ULIBS1D = $(DIR_BASE)/memdebug/d/libmemdebug.a \
  $(DIR_BASE)/graphics/1d/libgraphics.a $(DIR_BASE)/gui/d/libgui.a
ULIBS2D = $(DIR_BASE)/memdebug/d/libmemdebug.a \
  $(DIR_BASE)/graphics/2d/libgraphics.a $(DIR_BASE)/gui/d/libgui.a
GLIBS = `pkg-config --libs gtkglext-x11-1.0`

INCLUDE_OPTION = -I$(DIR_BASE)/memdebug -I$(DIR_BASE)/graphics \
  -I$(DIR_BASE)/gui `pkg-config --cflags gtkglext-x11-1.0`

HW20_FILES = CompressedSparsityPattern.o LaDouble.o Matrix.o NumPtr.o \
  Precondition.o SparseMatrix.o SparsityPattern.o Solver.o Vector.o \
  dhels.o dheqr.o dorth.o

all : iterative_improvement sor iterative_error_2d hw20

IterativeError2D.o :  IterativeError2D.C
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -DSPACEDIM=2 -o $*.o $<
%.o: %.C
	$(COMPILE.C++) $(INCLUDE_OPTION) -c -o $*.o $< 
%.o: %.f
	$(COMPILE.f) -c -o $*.o $<

iterative_improvement : IterativeImprovementMain.o iterative_improvement.o Level1D.o
	$(C++) -o $@ IterativeImprovementMain.o iterative_improvement.o \
	  Level1D.o $(ULIBS1D) $(GLIBS) -lblas

sor : SORMain.o iterative_improvement.o
	$(C++) -o $@ SORMain.o iterative_improvement.o $(ULIBS1D) \
	  $(GLIBS)

iterative_error_2d : IterativeError2D.o Level2D.o iterative_improvement_2d.o dicf.o fortio.o rcmp.o rbug.o
	$(C++) -o $@ IterativeError2D.o Level2D.o iterative_improvement_2d.o dicf.o fortio.o rcmp.o rbug.o $(ULIBS2D) $(GLIBS) -lblas -lgfortran

hw20 : hw20.o $(HW20_FILES)
	$(C++) -o $@ hw20.o $(HW20_FILES) $(ULIBS1D) $(GLIBS) -llapack -lblas

clean :
	-rm *.o core*
	-rm iterative_improvement sor iterative_error_2d hw20
depends.gnu :
	touch depends.gnu
	-$(RM) depends.gnu
	cat *.d > depends.gnu

include depends.gnu
