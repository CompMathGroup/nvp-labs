RM=rm -f
OPT=-O3
CFLAGS= $(OPT) -g -Wall -funroll-loops
CXXFLAGS=$(CFLAGS)
AR?=ar
TARGET=$(PREFIX)/libheatsolver.a
SRC= \
	heatsolver.cpp

OBJS= \
	$(PREFIX)/heatsolver.o

.PHONY: clean all

$(PREFIX)/%.o : %.cpp
	@echo "[CXX] "$<" => "$@
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "[AR] ... => "$(TARGET)
	@$(AR) rcs $@ $^

clean:
	$(RM) $(OBJS) $(TARGET)

$(PREFIX)/linalg.o:: linalg.cpp linalg.h array.h

$(PREFIX)/heatsolver.o:: heatsolver.cpp array.h coeff.h heatsolver.h solver.h
