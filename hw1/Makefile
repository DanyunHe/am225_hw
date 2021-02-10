# Include the configuration file from the parent directory
include ../config.mk

# Lists of executables
execs=casino cell prime sudoku 
all: $(execs)

# A Makefile target to remove all the built files
clean:
	rm -f $(execs)

# Makefile rule to build an executable. Note $@ is a shorthand for the file to
# be built and $^ is a shorthand for all the dependencies. (Can also use $< to
# mean the first dependency only.)
casino: casino.cc
	$(cxx) $(cflags) -o $@ $^

cell:cell.cc
	$(cxx) $(cflags) -o $@ $^

prime:prime.cc
	$(cxx) $(cflags) -o $@ $^

sudoku:sudoku.cc
	$(cxx) $(cflags) -o $@ $^

# The PHONY command tells GNU Make than "clean" is a phony target, which means
# that it doesn't correspond to an actual file called "clean'.
.PHONY: clean