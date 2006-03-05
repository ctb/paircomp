all: test-empty test-transit

test:
	./test-empty
	./test-transit

test-empty: test-empty.cc ../libpaircomp.a
	g++ -I../ -o test-empty test-empty.cc -L../ -lpaircomp

test-transit: test-transit.cc ../libpaircomp.a
	g++ -I../ -o test-transit test-transit.cc -L../ -lpaircomp
