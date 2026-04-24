CC=g++
CFLAGS=-c #-std=c++11
SOURCES= NC_3dFunctions.cpp  NC_CommonFunctions.cpp  NC_Input.cpp \
         NC_PostProcessing.cpp NC_Addresses.cpp    NC_EquationSystem.cpp \
         NC_Main.cpp NC_MLFMM.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../bin/NumCalc
GSLFLAG=-DUSE_GSL
LAPACKFLAG=-DUSE_LAPACK
OPTIM=-O2
GSL_PREFIX ?= /usr/local
LDFLAGS=-llapacke -lgsl -lblas -L$(GSL_PREFIX)/lib
CFLAGS+=-I$(GSL_PREFIX)/include
#DEBUG=-g

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) $(OPTIM) $(GSLFLAG) $(LAPACKFLAG) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(DEBUG) $(OPTIM) $(GSLFLAG) $(LAPACKFLAG) $<
clean:
	rm *.o
