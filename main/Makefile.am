bin_PROGRAMS = test test_ET simfit_example superjack_example
test_SOURCES = test.C
test_ET_SOURCES = test_ET.C
simfit_example_SOURCES = simfit_example.C
superjack_example_SOURCES = superjack_example.C

AM_CPPFLAGS = -I$(srcdir)/../include -fopenmp -g --std=c++14 @PYTHON_CFLAGS@

AM_LDFLAGS = -fopenmp -g 
LDADD = ../src/libcpsfit.a -lgsl -lgslcblas -lboost_serialization -lboost_python-py27 @PYTHON_LIBS@

