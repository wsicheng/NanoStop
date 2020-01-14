#
# stuff to make
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	build
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)'

ifndef VERBOSE
  QUIET := @
endif

#ROOFITINCLUDE =
#ifdef CMSSW_VERSION
#	ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
#endif

CC = g++
SNTCORE = ../NanoCORE
INCLUDE = $(shell root-config --cflags) -I$(SNTCORE)
CFLAGS = -Wall -Wno-unused-function -g -O3 -fPIC $(INCLUDE) $(EXTRACFLAGS)
ROOTLIBS = $(shell root-config --ldflags --cflags --libs) -lTMVA #-lEG -lGenVector
DICTINCLUDE = $(ROOTSYS)/include/Math/QuantFuncMathCore.h $(ROOTSYS)/include/TLorentzVector.h $(ROOTSYS)/include/Math/Vector4D.h

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs) -lEG -lGenVector -lTMVA -lRooFit

DIR = ./

SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
LIB = libStopLooper.so

COREDIR = ../NanoCORE
CORESOURCES = $(wildcard $(COREDIR)/*.cc)
COREOBJECTS = $(CORESOURCES:.cc=.o)
CORELIB = $(COREDIR)/NANO_CORE.so

# TOOLSSOURCES = ../CORE/Tools/goodrun.cc ../CORE/Tools/dorky/dorky.cc
# TOOLSOBJECTS = $(TOOLSSOURCES:.cc=.o)
# TOOLSLIB = libCORETools.so

LIBS = $(LIB) $(CORELIB)

EXE = runStopLooper

#
# how to make it
#

$(CORELIB): $(CORESOURCES)
	$(QUIET) echo "Linking $@"; \
	echo "Using Makefile from NanoCORE using $(MAKE) -j12"; \
	cd $(COREDIR) && $(MAKE) -j12 && cd $(DIR);

$(LIB):	$(OBJECTS)
	$(QUIET) echo "Linking $@"; \
	echo "$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)"; \
	$(LINKER) -shared -o $@ $(OBJECTS) $(DICT) $(LINKERFLAGS)

# the option "-Wl,-rpath,./" adds ./ to the runtime search path for libraries
$(EXE): $(LIBS)
	$(QUIET) echo "Building $@"; \
	echo "$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./"; \
	$(CC) -o $@ $(LIBS) $(ROOTLIBS) -Wl,-rpath,./

%.o: 	%.cc %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.C %.h
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.C
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o:    %.cxx
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

libs:	$(LIBS)

build:  $(EXE)

b: build

clean:
	rm -f \
	*.o \
	*.so \
	$(EXE)

.PHONY: clear cleanlibs

cleanlibs:
	rm -f \
	$(LIBS) \
	../NanoCORE/*.o \
	../NanoCORE/*.d \
	../NanoCORE/*.so \

clear: clean cleanlibs

endif
