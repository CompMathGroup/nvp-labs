OPT=-O3
CXXFLAGS=-g -Wall -I../core $(OPT) \
	-I../../include/$(ARCH)
MPREP=../../bin/$(ARCH)/mprep
LDFLAGS= \
	-L$(PREFIX) -L../../lib \

LIBS+=-lheatsolver

TARGET=$(PREFIX)/HeatLink.exe
OBJS=\
	$(PREFIX)/gate.o \
	$(PREFIX)/HeatLink.o

all: $(TARGET)

%.cpp : %.tm
	@echo "[MPREP] "$<
	@$(MPREP) $< > $@ || rm -f $@

$(PREFIX)/%.o : %.cpp
	@echo "[CXX] "$<" => "$@
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

$(TARGET): $(OBJS) $(PREFIX)/libheatsolver.a
	@echo "[LD] ... => "$@
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
