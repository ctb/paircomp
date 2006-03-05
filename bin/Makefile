INCLUDES=-I../lib/
CXXFLAGS=-g -Wall $(INCLUDES)

all: paircomp seqcomp find_patch window-extract

install: paircomp seqcomp find_patch
	cp paircomp seqcomp find_patch $(BINDIR)

paircomp: paircomp.cc ../lib/libpaircomp.a
	$(CXX) $(CXXFLAGS) -o paircomp paircomp.cc -L../lib -lpaircomp

seqcomp: seqcomp.cc ../lib/libpaircomp.a
	$(CXX) $(CXXFLAGS) -o seqcomp seqcomp.cc -L../lib -lpaircomp

find_patch: find_patch.cc ../lib/libpaircomp.a
	$(CXX) $(CXXFLAGS) -o find_patch find_patch.cc -L../lib -lpaircomp

window-extract: window-extract.cc ../lib/libpaircomp.a
	$(CXX) $(CXXFLAGS) -o window-extract window-extract.cc -L../lib -lpaircomp

clean:
	rm -f paircomp seqcomp find_patch window-extract
