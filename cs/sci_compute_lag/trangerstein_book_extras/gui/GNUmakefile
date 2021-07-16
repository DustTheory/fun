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

#INCLUDE_OPTION = -I$(DIR_BASE)/memdebug `pkg-config --cflags gtkglext-x11-1.0`
INCLUDE_OPTION = -I$(DIR_BASE)/memdebug -I$(DIR_BASE)/graphics `pkg-config --cflags gtkglext-x11-1.0`
ULIBS = $(DIR_BASE)/memdebug/d/libmemdebug.a $(DIR_BASE)/graphics/1d/libgraphics.a -lgthread-2.0
GLIBS = `pkg-config --libs gtkglext-x11-1.0`

DEBUGCOMPILEC = $(CPLUSPLUS) $(DEBUGFLAGS) -DDEBUG
DEBUGCOMPILEf = $(FC) $(DEBUGFFLAGS)
DEBUGCOMPILEc = $(CC) $(DEBUGFLAGS)
OPTCOMPILEC = $(CPLUSPLUS) $(OPTFLAGS) 
OPTCOMPILEf = $(FC) $(OPTFFLAGS)
OPTCOMPILEc = $(CC) $(OPTFLAGS)


CSOURCES = NISLList.C VirtualInput.C InputParameter.C \
  GTKGUIVirtualInput.C GUIInputParameter.C GUIEnumInputParameter.C \
  Thread.C GUIBase.C \
  GTKGUI.C GUI.C
fSOURCES =
cSOURCES =

CDEBUGOBJECTS = $(foreach SOURCE, $(CSOURCES), d/$(subst .C,.o,$(SOURCE)))
COPTOBJECTS = $(foreach SOURCE, $(CSOURCES), o/$(subst .C,.o,$(SOURCE)))

d/%.o: %.C
	$(DEBUGCOMPILEC) $(INCLUDE_OPTION) -c -o d/$*.o $<
d/%.o: %.f
	$(DEBUGCOMPILEf) -c -o d/$*.o $<
o/%.o: %.C
	$(OPTCOMPILEC) $(INCLUDE_OPTION) -c -o o/$*.o $<
o/%.o: %.f
	$(OPTCOMPILEf) -c -o o/$*.o $<

all : d/libgui.a o/libgui.a
	-rm depends.gnu
	cat d/*.d > depends.gnu
d/libgui.a : $(CDEBUGOBJECTS)
	ar -rv d/libgui.a $(CDEBUGOBJECTS)
o/libgui.a : $(COPTOBJECTS)
	ar -rv o/libgui.a $(COPTOBJECTS)

programs : testThread guiivp
testThread : d/testThread.o
	$(LINKC) -o $@ d/testThread.o d/libgui.a $(ULIBS) 
guiivp : d/guiivp.o d/guimain.o d/integrate.o
	$(LINKC) -o $@ d/guiivp.o d/guimain.o d/integrate.o \
	  d/libgui.a $(ULIBS) $(GLIBS)

clean :
	-rm d/*.o d/*.a d/*.d
	-rm o/*.o o/*.a o/*.d
	-rm guiivp testThread

zip :
	-tar -cvf tarfile README *.[CHchif] GNUmakefile
depends.gnu :
	touch depends.gnu
	-mkdir d
	-mkdir o
