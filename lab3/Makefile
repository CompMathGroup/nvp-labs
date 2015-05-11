BISON=bison
RAGEL=ragel
CXX=g++ -fmax-errors=5
#CXX=clang++ -ferror-limit=5
CXXFLAGS=-std=c++11 -g -O3 -Wall -I/usr/include/eigen3
OBJS=Compiler.o Polygon.o Simplex.o

%.tab.cc %.tab.hh : %.yy
	${BISON} -t -d --report=state -L c++ $<

%.lex.cc : %.rl
	${RAGEL} -G2 -o $@ $<

all: main

clean:
	rm -f *.tab.cc *.lex.cc *.hh *.o *.output main

parser: config.tab.o

lexer: config.lex.o

main: parser lexer main.o ${OBJS}
	${CXX} -o $@ config.tab.o config.lex.o main.o ${OBJS}

config.tab.o :: config.yy Lexer.h config.tab.hh Value.h Reg.h
config.lex.o :: config.rl config.yy Lexer.h config.tab.hh Value.h
config.tab.o :: config.yy Lexer.h config.tab.hh Value.h Reg.h

Compiler.o :: Value.h Reg.h
Polygon.o :: Value.h Reg.h
main.o :: Value.h Reg.h Container.h Simplex.h
Simlpex.o :: Simplex.h
