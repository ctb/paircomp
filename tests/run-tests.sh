#! /bin/bash
./test-python-interface.py && \
./test-regress.py && \
./test-transitive.py && \
./test-transitive-2.py && \
./test-compare-mussa.py && \
./test-errors.py
/bin/rm -f *.cmp
