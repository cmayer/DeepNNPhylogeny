## This is a bug fix version of the original quartet pattern counter. A more efficient version will be made avaiable soon.
## The executable is still called quartet-pattern-counter-v1.1.
# CFLAGS  = -g3 -ggdb -Wall
CFLAGS = -O3 # -Wall -Wextra

INCL =    -I .

SRC    = quartet-pattern-counter.cpp

HEADER = numpy.h symbol-cartesian-product_string.h



all:    quartet-pattern-counter-v1.1

quartet-pattern-counter-v1.1: $(SRC) $(HEADER)
	g++ $(CFLAGS) $(INCL) $(SRC) -o quartet-pattern-counter-v1.1


clean:
	rm -f quartet-pattern-counter-v1.1

