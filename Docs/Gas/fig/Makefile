.PHONY: all

SRC=gas.mp
RES=$(SRC:.mp=-0.eps)

MPOST=mpost

%-0.eps: %.mp
	$(MPOST) $<

all: $(RES)
