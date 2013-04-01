RM=rm -f
OPT=-O3
CFLAGS= $(OPT) -g -Wall -funroll-loops
CXXFLAGS=$(CFLAGS)
AR=ar -rcs
TARGET=$(PREFIX)/libgassolver.a
SRC= \
	gdsolver.cpp \
	advsolver.cpp \
	riemann.cpp

OBJS= \
	$(PREFIX)/gdsolver.o \
	$(PREFIX)/advsolver.o \
	$(PREFIX)/riemann.o

.PHONY: clean all

$(PREFIX)/%.o : %.cpp
	@echo "[CXX] "$<" => "$@
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(TARGET)

$(TARGET): $(OBJS)
	@echo "[AR] ... => "$(TARGET)
	@$(AR) -o $@ $^

clean: 
	$(RM) $(OBJS) $(TARGET)

linalg.o:: linalg.cpp linalg.h array.h 

solver.o:: solver.cpp solver.h array.h 
