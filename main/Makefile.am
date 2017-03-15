bin_PROGRAMS = test
test_SOURCES = test.C

AM_CPPFLAGS = -I$(srcdir)/../include -fopenmp -g --std=c++11 @PYTHON_CFLAGS@

AM_LDFLAGS = -fopenmp -g 
LDADD = ../src/libcpsfit.a -lgsl -lgslcblas -lboost_serialization -lboost_python-py27 @PYTHON_LIBS@

