RM=rm -f
OPT=-O3
CFLAGS= $(OPT) -g -Wall -funroll-loops
CXXFLAGS=$(CFLAGS)
AR=ar -rcs
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
	@$(AR) -o $@ $^

clean: 
	$(RM) $(OBJS) $(TARGET)
