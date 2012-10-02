CC=gcc
CXX=g++
RM=rm -f
CFLAGS=-O0 -g
CXXFLAGS=-O0 -g
LDFLAGS=-lm
TARGET=gas
OBJS=gas.o matrix.o

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

clean: 
	$(RM) $(OBJS) $(TARGET)

gas.o: matrix.h

matrix.o: matrix.h
