CC=gcc
CXX=clang++ -std=c++11
RM=rm -f
CFLAGS=-O0 -g -Wall
CXXFLAGS=$(CFLAGS)
LDFLAGS=-g -lm
TARGET=gas
OBJS=main.o gas.o solver.o

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

purge:
	$(RM) sod*.csv

clean: purge
	$(RM) $(OBJS) $(TARGET)

gas.o:: coeff.h 

main.o:: main.h array.h 

solver.o:: solver.h array.h 
