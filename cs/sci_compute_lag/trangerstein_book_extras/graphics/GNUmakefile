SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))

DEBUGFLAGS = -g
OPTFLAGS = -O3
CPLUSPLUS = g++ -Wnon-virtual-dtor -Woverloaded-virtual -g -MMD -DUSE_GTK -DNO_THREAD
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -Wall
CC = gcc -g -MMD
FC = gfortran
LINKC = $(CPLUSPLUS)
INCLUDE_OPTION = -I$(DIR_BASE)/memdebug `pkg-config --cflags gtkglext-x11-1.0`
#INCLUDE_OPTION = `pkg-config --cflags gtkglext-x11-1.0`

DEBUGCOMPILEC = $(CPLUSPLUS) $(DEBUGFLAGS) -DDEBUG
DEBUGCOMPILEf = $(FC) $(DEBUGFFLAGS)
DEBUGCOMPILEc = $(CC) $(DEBUGFLAGS)
OPTCOMPILEC = $(CPLUSPLUS) $(OPTFLAGS) 
OPTCOMPILEf = $(FC) $(OPTFFLAGS)
OPTCOMPILEc = $(CC) $(OPTFLAGS)

#DEBUGLIBS = $(DIR_BASE)/memdebug/d/libmemdebug.a -lgthread-2.0
#OPTLIBS = $(DIR_BASE)/memdebug/o/libmemdebug.a -lgthread-2.0
#DEBUGLIBS = -lgthread-2.0
#OPTLIBS = -lgthread-2.0
DEBUGLIBS = $(DIR_BASE)/memdebug/d/libmemdebug.a
OPTLIBS = $(DIR_BASE)/memdebug/o/libmemdebug.a
GLIBS = `pkg-config --libs gtkglext-x11-1.0`

GTK_SOURCES = GTKColormap.C GTKWindow.C GDKColorArray.C GTKGLWindow.C
CSOURCES = ISLList.C NumPtr.C \
  Palette.C VGT.C XYGraphTool.C AreaGraphTool.C \
  VirtualGLWindow.C VolGraphTool.C $(GTK_SOURCES)
cSOURCES = 
fSOURCES =

CDEBUGOBJECTS1 = $(foreach SOURCE, $(CSOURCES), 1d/$(subst .C,.o,$(SOURCE)))
CDEBUGOBJECTS2 = $(foreach SOURCE, $(CSOURCES), 2d/$(subst .C,.o,$(SOURCE)))
CDEBUGOBJECTS3 = $(foreach SOURCE, $(CSOURCES), 3d/$(subst .C,.o,$(SOURCE)))
COPTOBJECTS1 = $(foreach SOURCE, $(CSOURCES), 1o/$(subst .C,.o,$(SOURCE)))
COPTOBJECTS2 = $(foreach SOURCE, $(CSOURCES), 2o/$(subst .C,.o,$(SOURCE)))
COPTOBJECTS3 = $(foreach SOURCE, $(CSOURCES), 3o/$(subst .C,.o,$(SOURCE)))

1d/%.o: %.C
	$(DEBUGCOMPILEC) -DSPACEDIM=1 $(INCLUDE_OPTION) -c -o 1d/$*.o $< 
2d/%.o: %.C
	$(DEBUGCOMPILEC) -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o 2d/$*.o $< 
3d/%.o: %.C
	$(DEBUGCOMPILEC) -DSPACEDIM=3 $(INCLUDE_OPTION) -c -o 3d/$*.o $< 
1o/%.o: %.C
	$(OPTCOMPILEC) -DSPACEDIM=1 $(INCLUDE_OPTION) -c -o 1o/$*.o $< 
2o/%.o: %.C
	$(OPTCOMPILEC) -DSPACEDIM=2 $(INCLUDE_OPTION) -c -o 2o/$*.o $< 
3o/%.o: %.C
	$(OPTCOMPILEC) -DSPACEDIM=3 $(INCLUDE_OPTION) -c -o 3o/$*.o $< 

all : 1d/libgraphics.a 2d/libgraphics.a 3d/libgraphics.a \
  1o/libgraphics.a 2o/libgraphics.a 3o/libgraphics.a
	-rm depends.gnu
	-cat 2d/*.d 2o/*.d > depends.gnu
1d/libgraphics.a : $(CDEBUGOBJECTS1)
	ar -rv 1d/libgraphics.a $(CDEBUGOBJECTS1)
2d/libgraphics.a : $(CDEBUGOBJECTS2)
	ar -rv 2d/libgraphics.a $(CDEBUGOBJECTS2)
3d/libgraphics.a : $(CDEBUGOBJECTS3)
	ar -rv 3d/libgraphics.a $(CDEBUGOBJECTS3)
1o/libgraphics.a : $(COPTOBJECTS1)
	ar -rv 1o/libgraphics.a $(COPTOBJECTS1)
2o/libgraphics.a : $(COPTOBJECTS2)
	ar -rv 2o/libgraphics.a $(COPTOBJECTS2)
3o/libgraphics.a : $(COPTOBJECTS3)
	ar -rv 3o/libgraphics.a $(COPTOBJECTS3)

programs : testXYGraphTool newton gasDynamicsRiemannProblem \
  testAreaGraphToolColor testAreaGraphToolContour testVolGraphTool2 \
  testVolGraphTool3 testVolGraphTool4 testVolGraphTool5 
testXYGraphTool : 1d/testXYGraphTool.o
	$(LINKC) -o $@ 1d/testXYGraphTool.o 1d/libgraphics.a $(DEBUGLIBS) \
	  $(GLIBS)
newton : 1d/newton.o
	$(LINKC) -o $@ 1d/newton.o 1d/libgraphics.a $(DEBUGLIBS) $(GLIBS)
gasDynamicsRiemannProblem : 1d/gasDynamicsRiemannProblem.o
	$(LINKC) -o $@ 1d/gasDynamicsRiemannProblem.o 1d/libgraphics.a \
	  $(DEBUGLIBS) $(GLIBS)
testAreaGraphToolColor : 2d/testAreaGraphToolColor.o
	$(LINKC) -o $@ 2d/testAreaGraphToolColor.o 2d/libgraphics.a \
	  $(DEBUGLIBS) $(GLIBS)
testAreaGraphToolContour : 2d/testAreaGraphToolContour.o
	$(LINKC) -o $@ 2d/testAreaGraphToolContour.o 2d/libgraphics.a \
	  $(DEBUGLIBS) $(GLIBS)
testVolGraphTool2 : 2d/testVolGraphTool2.o
	$(LINKC) -o $@ 2d/testVolGraphTool2.o 2d/libgraphics.a \
	  $(DEBUGLIBS) $(GLIBS)
testVolGraphTool3 : 3d/testVolGraphTool3.o
	$(LINKC) -o $@ 3d/testVolGraphTool3.o 3d/libgraphics.a \
	  $(DEBUGLIBS) $(GLIBS)
testVolGraphTool4 : 3d/testVolGraphTool4.o 3d/SurfaceDataSet.o
	$(LINKC) -o $@ 3d/testVolGraphTool4.o 3d/SurfaceDataSet.o \
	  3d/libgraphics.a $(DEBUGLIBS) $(GLIBS)
testVolGraphTool5 : 3d/testVolGraphTool5.o 3d/SurfaceDataSet.o
	$(LINKC) -o $@ 3d/testVolGraphTool5.o 3d/SurfaceDataSet.o \
	  3d/libgraphics.a $(DEBUGLIBS) $(GLIBS)

clean :
	-rm 1d/*.o 1d/*.a 1d/*.d
	-rm 2d/*.o 2d/*.a 2d/*.d
	-rm 3d/*.o 3d/*.a 3d/*.d
	-rm 1o/*.o 1o/*.a 1o/*.d
	-rm 2o/*.o 2o/*.a 2o/*.d
	-rm 3o/*.o 3o/*.a 3o/*.d
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
