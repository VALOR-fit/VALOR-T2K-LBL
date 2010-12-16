ROOTLIBS   = $(shell root-config --glibs) #-lNew
ROOTCFLAGS = $(shell root-config --cflags)

CXXFLAGS   = -O3 -g -fPIC $(ROOTCFLAGS)
SYSLIB     = -lm

LINK_ARGS_BIN = $(SYSLIB) $(ROOTLIBS)

BIN = ../bin/test
OBJS = test.o HadErrCalc.o

all: test

test: Makefile $(OBJS)
	g++ $(CXXFLAGS)  -o $(BIN)  $(OBJS) $(LINK_ARGS_BIN)

.cc.o:
	g++ $(CXXFLAGS)  -c $<

test.o: HadErrCalc.h
HadErrCalc.o: HadErrCalc.h

clean:
	rm $(PROGS) *~ *.o $(BIN)

