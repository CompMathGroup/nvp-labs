CC=gcc 
CXX=g++ 
RM=rm -f
CFLAGS=-m32 -O3 -g -Wall
CXXFLAGS=$(CFLAGS)
LDFLAGS=-g -lm -lprofiler
AR=ar -rcs
TARGET=libgassolver32.a
OBJS=solver.o linalg.o riemann.o

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) -o $@ $^
#	$(CXX) $(LDFLAGS) -o $@ $^

purge:
	$(RM) sod*.csv

clean: purge
	$(RM) $(OBJS) $(TARGET)

linalg.o:: linalg.cpp linalg.h array.h 

solver.o:: solver.cpp solver.h array.h 
