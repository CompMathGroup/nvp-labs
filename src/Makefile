CC=gcc
CXX=clang++ 
#-std=c++11
RM=rm -f
CFLAGS=-O3 -g -Wall
CXXFLAGS=$(CFLAGS)
LDFLAGS=-g -lm -lprofiler
TARGET=gas
OBJS=main.o gas.o solver.o schemes.o

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

purge:
	$(RM) sod*.csv

clean: purge
	$(RM) $(OBJS) $(TARGET)

gas.o:: main.cpp coeff.h 

main.o:: main.cpp main.h array.h 

schemes.o:: schemes.cpp schemes.h coeff.h 

solver.o:: solver.cpp solver.h array.h 
