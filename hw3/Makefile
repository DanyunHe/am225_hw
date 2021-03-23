# Include the configuration file
include ../config.mk

# Lists of files to be built
objs=file_output.o conj_grad.o rbf.o lp_solve.o
src=$(patsubst %.o,%.cc,$(objs))
execs=rbf_time matrix_mul

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

liblinal.a: $(objs)
	rm -f $@
	ar rs $@ $^

%.o: %.cc
	$(cxx) $(cflags) -c $<

rbf_time: rbf_time.cc liblinal.a
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

matrix_mul: matrix_mul.cc
	$(cxx) $(cflags) -o $@ $^ $(lp_lflags)

.PHONY: clean depend