SD = $(shell pwd)
DIR_BASE = $(shell dirname $(SD))

DEBUGFLAGS = -g
OPTFLAGS = -O3
CPLUSPLUS = g++ -Wnon-virtual-dtor -Woverloaded-virtual -g -MMD
DEBUGFFLAGS = -g -u -pedantic-errors -fsecond-underscore -fimplicit-none -Wsurprising -Werror
OPTFFLAGS = -O3 -g -Wall
CC = gcc -g -MMD
FC = gfortran
LINKC = $(CPLUSPLUS)

DEBUGCOMPILEC = $(CPLUSPLUS) $(DEBUGFLAGS) -DDEBUG
DEBUGCOMPILEf = $(FC) $(DEBUGFFLAGS)
DEBUGCOMPILEc = $(CC) $(DEBUGFLAGS)
OPTCOMPILEC = $(CPLUSPLUS) $(OPTFLAGS)
OPTCOMPILEf = $(FC) $(OPTFFLAGS)
OPTCOMPILEc = $(CC) $(OPTFLAGS)

CSOURCES = Debug.C MemoryDebugger.C Tracer.C TimedObject.C
cSOURCES = second.c
fSOURCES = testtypes.f

CDEBUGOBJECTS = $(foreach SOURCE, $(CSOURCES), d/$(subst .C,.o,$(SOURCE)))
fDEBUGOBJECTS = $(foreach SOURCE, $(fSOURCES), d/$(subst .f,.o,$(SOURCE)))
cDEBUGOBJECTS = $(foreach SOURCE, $(cSOURCES), d/$(subst .c,.o,$(SOURCE)))
COPTOBJECTS = $(foreach SOURCE, $(CSOURCES), o/$(subst .C,.o,$(SOURCE)))
fOPTOBJECTS = $(foreach SOURCE, $(fSOURCES), o/$(subst .f,.o,$(SOURCE)))
cOPTOBJECTS = $(foreach SOURCE, $(cSOURCES), o/$(subst .c,.o,$(SOURCE)))

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

all : d/libmemdebug.a o/libmemdebug.a
	-rm depends.gnu
	-cat d/*.d o/*.d > depends.gnu
d/libmemdebug.a : $(CDEBUGOBJECTS) $(fDEBUGOBJECTS) $(cDEBUGOBJECTS)
	ar -rv d/libmemdebug.a $(CDEBUGOBJECTS) $(fDEBUGOBJECTS) \
	  $(cDEBUGOBJECTS)
o/libmemdebug.a : $(COPTOBJECTS) $(fOPTOBJECTS) $(cOPTOBJECTS)
	ar -rv o/libmemdebug.a $(COPTOBJECTS) $(fOPTOBJECTS) \
	  $(cOPTOBJECTS)
programs : testMemoryDebugger testTimedObject
testMemoryDebugger : d/testMemoryDebugger.o d/testfortroutine.o \
  d/libmemdebug.a
	$(LINKC) -o $@ d/testMemoryDebugger.o d/testfortroutine.o \
	  d/libmemdebug.a
testTimedObject : d/testTimedObject.o d/libmemdebug.a
	$(LINKC) -o $@ d/testTimedObject.o d/libmemdebug.a

clean :
	-rm d/*.o d/*.a d/*.d
	-rm o/*.o o/*.a o/*.d
zip :
	-tar -cvf tarfile *.[CHchif] GNUmakefile
depends.gnu :
	touch depends.gnu
	-mkdir d
	-mkdir o

include depends.gnu
