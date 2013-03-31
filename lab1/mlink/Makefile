CXX=g++
OPT=-O3
CXXFLAGS= -I../core $(OPT) -g \
	-I/usr/local/Wolfram/Mathematica/8.0/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/

LDFLAGS=-static-libgcc -L../core
MCC=/usr/local/Wolfram/Mathematica/8.0/Executables/mcc

all: unix32 unix64

gate32.o : gate.cpp
	$(CXX) -m32 $(CXXFLAGS) -c -o $@ $<

gate64.o : gate.cpp
	$(CXX) -m64 $(CXXFLAGS) -c -o $@ $<

unix32: GasLink.tm gate32.o
	mkdir -p Unix_x86
	$(MCC) $^ -b32 -o Unix_x86/GasLink.exe $(LDFLAGS) -lgassolver32
	
unix64: GasLink.tm gate64.o
	mkdir -p Unix_x86-64
	$(MCC) $^ -b64 -o Unix_x86-64/GasLink.exe $(LDFLAGS) -lgassolver64

clean:
	rm -f *.o
