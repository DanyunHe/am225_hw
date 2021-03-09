# Include the configuration file from the parent directory
include ../config.mk

# Lists of executables
objs=combinations.o fsal.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs=tree fsal_test
all: $(execs)

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(execs)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<
# Makefile rule to build an executable. Note $@ is a shorthand for the file to
# be built and $^ is a shorthand for all the dependencies. (Can also use $< to
# mean the first dependency only.)
tree: tree.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

fsal_test: fsal_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

# combinations: combinations.cc
# 	$(cxx) $(cflags) -o $@ $^
# The PHONY command tells GNU Make than "clean" is a phony target, which means
# that it doesn't correspond to an actual file called "clean'.
.PHONY: clean