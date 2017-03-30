bin_PROGRAMS = test test_ET
test_SOURCES = test.C
test_ET_SOURCES = test_ET.C

AM_CPPFLAGS = -I$(srcdir)/../include -fopenmp -g --std=c++14 @PYTHON_CFLAGS@

AM_LDFLAGS = -fopenmp -g 
LDADD = ../src/libcpsfit.a -lgsl -lgslcblas -lboost_serialization -lboost_python-py27 @PYTHON_LIBS@

