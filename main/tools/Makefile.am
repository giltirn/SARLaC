exampledir = $(prefix)/example
testdir = $(prefix)/test
bin_PROGRAMS = gaussian_jackknife hdf5_calc hdf5_combine hdf5_print hdf5_superjackknife_boost
gaussian_jackknife_SOURCES = gaussian_jackknife.C
hdf5_calc_SOURCES = hdf5_calc.C
hdf5_combine_SOURCES = hdf5_combine.C
hdf5_print_SOURCES = hdf5_print.C
hdf5_superjackknife_boost_SOURCES = hdf5_superjackknife_boost.C
AM_CPPFLAGS = -I$(top_srcdir)/include -I$(srcdir)/include
LDADD =$(top_builddir)/src/libsarlac.a
