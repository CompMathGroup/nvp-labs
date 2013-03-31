CC=gcc 
CXX=g++ 
RM=rm -f
OPT=-O3
CFLAGS32=-m32 $(OPT) -g -Wall -funroll-loops
CXXFLAGS32=$(CFLAGS32)
CFLAGS64=-m64 $(OPT) -g -Wall -funroll-loops
CXXFLAGS64=$(CFLAGS64)
LDFLAGS=-g -lm -lprofiler
AR=ar -rcs
TARGET32=libgassolver32.a
TARGET64=libgassolver64.a
SRC=gdsolver.cpp advsolver.cpp riemann.cpp
OBJS32=$(SRC:.cpp=32.o)
OBJS64=$(SRC:.cpp=64.o)

.PHONY: clean all

%32.o : %.cpp
	$(CXX) $(CXXFLAGS32) -c -o $@ $<

%64.o : %.cpp
	$(CXX) $(CXXFLAGS64) -c -o $@ $<

all: $(TARGET32) $(TARGET64)

$(TARGET32): $(OBJS32)
	$(AR) -o $@ $^

$(TARGET64): $(OBJS64)
	$(AR) -o $@ $^

clean: 
	$(RM) $(OBJS32) $(OBJS64) $(TARGET32) $(TARGET64)

linalg.o:: linalg.cpp linalg.h array.h 

solver.o:: solver.cpp solver.h array.h 
