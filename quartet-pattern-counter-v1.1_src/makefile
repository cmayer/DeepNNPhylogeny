# CFLAGS  = -g3 -ggdb -Wall
CFLAGS = -O3 # -Wall -Wextra

INCL =    -I .

SRC    = quartet-pattern-counter.cpp symbol-cartesian-product_string.h

HEADER = numpy.h



all:    quartet-pattern-counter-v1.1

quartet-pattern-counter-v1.1: $(SRC) $(HEADER)
	g++ $(CFLAGS) $(INCL) $(SRC) -o quartet-pattern-counter-v1.1


clean:
	rm -f quartet-pattern-counter-v1.1

