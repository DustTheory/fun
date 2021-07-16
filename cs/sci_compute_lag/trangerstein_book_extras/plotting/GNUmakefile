SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))

DEBUGFLAGS = -g -DNO_COLOR_ARRAY -DMGL_SRC -DNDEBUG -DMGL_STATIC_DEFINE
OPTFLAGS = -O3 -DNO_COLOR_ARRAY -DMGL_SRC -DNDEBUG -DMGL_STATIC_DEFINE
CPLUSPLUS = g++ -Wnon-virtual-dtor -Woverloaded-virtual -g -MMD -DUSE_GTK -DNO_THREAD
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -Wall
CC = gcc -g -MMD
FC = gfortran
LINKC = $(CPLUSPLUS)
INCLUDE_OPTION = -I$(DIR_BASE)/memdebug -I$(DIR_BASE)/mathgl-2.1.3.1/include -I$(DIR_BASE)/mathgl-2.1.3.1/include/mgl2

DEBUGCOMPILEC = $(CPLUSPLUS) $(DEBUGFLAGS) -DDEBUG
DEBUGCOMPILEf = $(FC) $(DEBUGFFLAGS)
DEBUGCOMPILEc = $(CC) $(DEBUGFLAGS)
OPTCOMPILEC = $(CPLUSPLUS) $(OPTFLAGS) 
OPTCOMPILEf = $(FC) $(OPTFFLAGS)
OPTCOMPILEc = $(CC) $(OPTFLAGS)

DEBUGLIBS = $(DIR_BASE)/memdebug/d/libmemdebug.a -lgthread-2.0
OPTLIBS = $(DIR_BASE)/memdebug/o/libmemdebug.a -lgthread-2.0
GLIBS = $(DIR_BASE)/mathgl-2.1.3.1/src/libmgl.so.7.0.0

%.o: %.C
	$(DEBUGCOMPILEC) $(INCLUDE_OPTION) -c -o $*.o $< 

all : testMGLContour testMGLColor testMGLSurface testMGLIsoSurface
	-rm depends.gnu
	-cat *.d > depends.gnu
#testXYGraphTool : 1d/testXYGraphTool.o
#	$(LINKC) -o $@ 1d/testXYGraphTool.o 1d/libgraphics.a $(DEBUGLIBS) \
#	  $(GLIBS)
#newton : 1d/newton.o
#	$(LINKC) -o $@ 1d/newton.o 1d/libgraphics.a $(DEBUGLIBS) $(GLIBS)
#gasDynamicsRiemannProblem : 1d/gasDynamicsRiemannProblem.o
#	$(LINKC) -o $@ 1d/gasDynamicsRiemannProblem.o 1d/libgraphics.a \
#	  $(DEBUGLIBS) $(GLIBS)
#testAreaGraphToolColor : 2d/testAreaGraphToolColor.o
#	$(LINKC) -o $@ 2d/testAreaGraphToolColor.o 2d/libgraphics.a \
#	  $(DEBUGLIBS) $(GLIBS)
testMGLContour : testMGLContour.o
	$(LINKC) -o $@ testMGLContour.o $(DEBUGLIBS) $(GLIBS)
testMGLColor : testMGLColor.o
	$(LINKC) -o $@ testMGLColor.o $(DEBUGLIBS) $(GLIBS)
testMGLSurface : testMGLSurface.o
	$(LINKC) -o $@ testMGLSurface.o $(DEBUGLIBS) $(GLIBS)
testMGLIsoSurface : testMGLIsoSurface.o
	$(LINKC) -o $@ testMGLIsoSurface.o $(DEBUGLIBS) $(GLIBS)
#testVolGraphTool2 : 2d/testVolGraphTool2.o
#	$(LINKC) -o $@ 2d/testVolGraphTool2.o 2d/libgraphics.a \
#	  $(DEBUGLIBS) $(GLIBS)
#testVolGraphTool3 : 3d/testVolGraphTool3.o
#	$(LINKC) -o $@ 3d/testVolGraphTool3.o 3d/libgraphics.a \
#	  $(DEBUGLIBS) $(GLIBS)
#testVolGraphTool4 : 3d/testVolGraphTool4.o 3d/SurfaceDataSet.o
#	$(LINKC) -o $@ 3d/testVolGraphTool4.o 3d/SurfaceDataSet.o \
#	  3d/libgraphics.a $(DEBUGLIBS) $(GLIBS)
#testVolGraphTool5 : 3d/testVolGraphTool5.o 3d/SurfaceDataSet.o
#	$(LINKC) -o $@ 3d/testVolGraphTool5.o 3d/SurfaceDataSet.o \
#	  3d/libgraphics.a $(DEBUGLIBS) $(GLIBS)

clean :
	-rm 1d/*.o 1d/*.a 1d/*.d
	-rm 2d/*.o 2d/*.a 2d/*.d
	-rm 3d/*.o 3d/*.a 3d/*.d
	-rm 1o/*.o 1o/*.a 1o/*.d
	-rm 2o/*.o 2o/*.a 2o/*.d
	-rm 3o/*.o 3o/*.a 3o/*.d
	-rm testMGLContour testMGLColor testMGLSurface testMGLIsoSurface
zip :
	-tar -cvf tarfile *.[CH] GNUmakefile
depends.gnu :
	touch depends.gnu
	-mkdir 1d
	-mkdir 2d
	-mkdir 3d
	-mkdir 1o
	-mkdir 2o
	-mkdir 3o

include depends.gnu
