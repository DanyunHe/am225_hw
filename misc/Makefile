# Load the common configuration file
include ../config.mk

objs=write_png.o
src=$(patsubst %.o,%.cc,$(objs))
execs=wpng_example paper_png

all: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

wpng_example: wpng_example.cc write_png.o
	$(cxx) $(cflags) $(iflags) $(png_iflags) -o $@ $^ $(lflags) $(png_lflags)

paper_png: paper_png.cc write_png.o
	$(cxx) $(cflags) $(iflags) $(png_iflags) -o $@ $^ $(lflags) $(png_lflags)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) $(png_iflags) -c $<

clean:
	rm -f $(execs) $(objs) libf2d.a

clean-all:
	$(MAKE) -C ../tgmg clean
	rm -f $(execs) $(objs) libf2d.a

.PHONY: clean all executables depend
