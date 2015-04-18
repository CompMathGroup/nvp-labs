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

$(PREFIX)/linalg.o:: linalg.cpp linalg.h array.h

$(PREFIX)/gdsolver.o:: gdsolver.cpp solver.h array.h coeff.h

$(PREFIX)/advsolver.o:: advsolver.cpp solver.h array.h coeff.h
