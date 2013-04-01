OPT=-O3
CXXFLAGS=-g -Wall -I../core $(OPT) \
	-I/usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/
MPREP=/usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/mprep
LDFLAGS= \
	-L$(PREFIX) -L../lib \
	-static-libstdc++ -static-libgcc

LIBS+=-lgassolver

TARGET=$(PREFIX)/GasLink.exe
OBJS=\
	$(PREFIX)/gate.o \
	$(PREFIX)/GasLink.o

all: $(TARGET)

%.cpp : %.tm
	@echo "[MPREP] "$<
	@$(MPREP) $< > $@

$(PREFIX)/%.o : %.cpp
	@echo "[CXX] "$<" => "$@
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

$(TARGET): $(OBJS)
	@echo "[LD] ... => "$@
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
