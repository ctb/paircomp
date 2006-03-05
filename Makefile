# system directory to put binaries into
BINDIR=/usr/local/bin

# Python command to use
PYTHON=python

### shouldn't need to change below this ###

all: lib-dir bin-dir python-dir

install: all
	cd bin/ && $(MAKE) BINDIR=$(BINDIR) install
	cd ../
	cd python/ && $(PYTHON) setup.py install
	cd ../

clean:
	cd lib/ && $(MAKE) clean
	cd bin/ && $(MAKE) clean
	cd python/ && rm -fr build
	cd ..

rebuild: python-dir-rebuild all

depend:
	cd lib/ && $(MAKE) depend
	cd ../

do-tests:
	cd tests && ./run-tests.sh

lib-dir:
	cd lib/ && $(MAKE) all
	cd ../

bin-dir:
	cd bin/ && $(MAKE) all
	cd ../

python-dir:
	cd python/ && $(PYTHON) setup.py build
	cd ../

python-dir-rebuild:
	cd python/ && $(PYTHON) setup.py clean --all
	cd ../
