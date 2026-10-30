CC ?= cc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra
LDLIBS = -lpthread -lm

# optional: make NATIVE=1 adds -march=native (faster bitsliced engine)
ifeq ($(NATIVE),1)
CFLAGS += -march=native
endif

BINS = bin/dltool bin/onestep bin/lintool bin/longdiff bin/lowgraph bin/lowgraph3 bin/hfcheck bin/krec bin/krsim bin/carrycheck bin/toysearch bin/toygraph bin/toymeasure

all: $(BINS)

bin/%: src/%.c src/dipper.h src/bitslice.h src/halffast.h src/graphalg.h src/bsinv.h src/krank.h src/toy32.h
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: all
	./bin/dltool test
	python3 steps/01_reference_check.py

clean:
	rm -rf bin

.PHONY: all test clean
