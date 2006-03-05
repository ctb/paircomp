CXX=g++
CXXFLAGS=-g -Wall

###### top level targets

all: libpaircomp.a

clean:
	rm -f *.o *.a

####### individual binary targets

COMMON_OBJS=MutableComparison.o ImmutableComparison.o \
	algorithms.o algorithms2.o

libpaircomp.a: $(COMMON_OBJS)
	ar rcv libpaircomp.a $(COMMON_OBJS)
	ranlib libpaircomp.a

####### various individual object targets

depend:
	makedepend -fMakefile.deps -- -I/usr/local/include $(CXXFLAGS) -- *.cc

include Makefile.deps
