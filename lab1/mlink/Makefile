OPT=-O3
CXXFLAGS=-s -Wall -I../core $(OPT) \
	-I../../include/$(ARCH)
MPREP=../../bin/$(ARCH)/mprep
LDFLAGS= \
	-L$(PREFIX) -L../../lib \

#-Wl,-dynamic-linker \

LIBS+=-lgassolver

TARGET=$(PREFIX)/GasLink.exe
OBJS=\
	$(PREFIX)/gate.o \
	$(PREFIX)/GasLink.o

all: $(TARGET)

%.cpp : %.tm
	@echo "[MPREP] "$<
	@$(MPREP) $< > $@ || rm -f $@

$(PREFIX)/%.o : %.cpp
	@echo "[CXX] "$<" => "$@
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

$(TARGET): $(OBJS) $(PREFIX)/libgassolver.a
	@echo "[LD] "$^" => "$@
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
