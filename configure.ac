#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.69])
AC_INIT([SARLaC],[1.0],[ckelly.physics@gmail.com])
AC_CONFIG_SRCDIR([include/utils.h])
AC_CONFIG_HEADERS([config.h])

CPPFLAGS="-Wno-deprecated-declarations $CPPFLAGS"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_REQUIRE_CPP

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi


# Check for c++14
CXXFLAGS="-std=c++14 $CXXFLAGS"
AC_MSG_CHECKING([for c++14 compatibility])
AC_COMPILE_IFELSE( [AC_LANG_PROGRAM([[]],
                                    [[]])],
                   [AC_MSG_RESULT([yes])],
                   [AC_MSG_FAILURE(["Could not compile with CXXFLAGS $CXXFLAGS"])]
                  )

# Check for openmp
AC_OPENMP
ac_openmp=no
if test "${OPENMP_CXXFLAGS}X" != "X"; then
  ac_openmp=yes
  CXXFLAGS="$OPENMP_CXXFLAGS $CXXFLAGS"
  LDFLAGS="$OPENMP_CXXFLAGS $LDFLAGS"
else
  AC_MSG_FAILURE(["OpenMP support is currently required"])
fi

# Checks for libraries.
PKG_PROG_PKG_CONFIG

#Check for python 2.7+
AC_ARG_WITH([python], AS_HELP_STRING([--with-python], [Build with Python]),[],[with_python=no])
if test "x$with_python" = "xyes"; then
   PKG_CHECK_MODULES([PYTHON], [python >= 2.7, python < 3.0 ], [AC_DEFINE([HAVE_PYTHON], [1], [Use Python])   ], [AC_MSG_FAILURE([Python version must be >=2.7 and <3.0])])
   LIBS="$PYTHON_LIBS $LIBS"
   CPPFLAGS="$PYTHON_CFLAGS $CPPFLAGS"
fi

m4_include([ax_boost_base.m4])
AX_BOOST_BASE([1.61],[], [AC_MSG_FAILURE([Require Boost version >=1.61])] )

CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS"
LDFLAGS="$BOOST_LDFLAGS $LDFLAGS"

#Check for boost
#AC_ARG_WITH([boost], AS_HELP_STRING([--with-boost], [Build with Boost]),[],[with_boost=no])
#if test "x$with_boost" = "xno"; then
#  AC_MSG_FAILURE("Must compile with boost. Use --with-boost=<install dir>", -1)
#elif test "x$with_boost" != "xyes"; then
#  CPPFLAGS="-I${with_boost}/include $CPPFLAGS"
#  LDFLAGS="-L${with_boost}/lib $LDFLAGS"
#fi

#AC_SUBST(BOOST_INSTALL_DIR,$boost_dir)

#Check for boost python
if test "x$with_python" = "xyes"; then
AC_CHECK_LIB([boost_python], [main], [LIBS="-lboost_python $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_python])])
fi

#Check for boost serialization
AC_CHECK_LIB([boost_serialization], [main], [LIBS="-lboost_serialization $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_serialization])])

#Check for boost timer
AC_CHECK_LIB([boost_timer], [main], [LIBS="-lboost_timer $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_timer])])

#Check for boost iostreams
AC_CHECK_LIB([boost_iostreams], [main], [LIBS="-lboost_iostreams $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_iostreams])])

#Check for boost filesystem
AC_CHECK_LIB([boost_filesystem], [main], [LIBS="-lboost_filesystem $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_filesystem])])

#Check for boost system
AC_CHECK_LIB([boost_system], [main], [LIBS="-lboost_system $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_system])])

#AM_CONDITIONAL([USE_PYTHON], [test "x$with_python" = "xyes"])

#Check for GSL
AC_ARG_WITH([gsl], AS_HELP_STRING([--with-gsl], [Build with GSL]),[],[with_gsl=no])
if test "x$with_gsl" = "xno"; then
  AC_MSG_FAILURE([GSL is required])
elif test "x$with_gsl" != "xyes"; then
  CPPFLAGS="-I${with_gsl}/include $CPPFLAGS"
  LDFLAGS="-L${with_gsl}/lib $LDFLAGS"
fi

AC_CHECK_LIB([gslcblas],[main], [LIBS="-lgslcblas $LIBS"], [AC_MSG_FAILURE([Could not link to gslcblas])])
AC_CHECK_LIB([gsl],[main], [LIBS="-lgsl $LIBS"], [AC_MSG_FAILURE([Could not link to gsl])])


#Check for HDF5
AC_ARG_WITH([hdf5], AS_HELP_STRING([--with-hdf5], [Build with HDF5]),[],[with_hdf5=no])
if test "x$with_hdf5" = "xyes"; then
   PKG_CHECK_MODULES([HDF5], [hdf5 >= 0 ], [AC_DEFINE([HAVE_HDF5], [1], [Use HDF5])   ], [AC_MSG_FAILURE([HDF5 version required])])
   LIBS="$HDF5_LIBS -lhdf5_cpp $LIBS"
   CPPFLAGS="$HDF5_CFLAGS $CPPFLAGS"
elif test "x$with_hdf5" != "xno"; then
  CPPFLAGS="-I${with_hdf5}/include $CPPFLAGS"
  LDFLAGS="-L${with_hdf5}/lib $LDFLAGS"
  AC_CHECK_LIB([hdf5],[main], [LIBS="-lhdf5 $LIBS"], [AC_MSG_FAILURE([Could not link to HDF5])])	
  AC_CHECK_LIB([hdf5_cpp],[main], [LIBS="-lhdf5_cpp $LIBS"], [AC_MSG_FAILURE([Could not link to HDF5 C++ library])])	
  AC_DEFINE([HAVE_HDF5], [1], [Use HDF5])
fi


#Check for Minuit2
AC_ARG_WITH([minuit2], AS_HELP_STRING([--with-minuit2], [Build with Minuit2]),[],[with_minuit2=no])
if test "x$with_minuit2" = "xyes"; then
   PKG_CHECK_MODULES([MINUIT2], [Minuit2 >= 0 ], [AC_DEFINE([HAVE_MINUIT2], [1], [Use Minuit2])   ], [AC_MSG_FAILURE([Minuit2 version required])])
   LIBS="$MINUIT2_LIBS -lMinuit2 $LIBS"
   CPPFLAGS="$MINUIT2_CFLAGS $CPPFLAGS"
elif test "x$with_minuit2" != "xno"; then
  CPPFLAGS="-I${with_minuit2}/include $CPPFLAGS"
  LDFLAGS="-L${with_minuit2}/lib $LDFLAGS"
  AC_CHECK_LIB([Minuit2],[main], [LIBS="-lMinuit2 $LIBS"], [AC_MSG_FAILURE([Could not link to Minuit2])])	
  AC_DEFINE([HAVE_MINUIT2], [1], [Use Minuit2])
fi



# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS(boost/format.hpp, [], [AC_MSG_FAILURE([["Could not compile with boost/format.hpp"]])] )
AC_CHECK_HEADERS(boost/spirit/home/x3.hpp,  [], [AC_MSG_FAILURE([["Could not compile with boost/spirit/home/x3.hpp"]])])

AC_PROG_RANLIB

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

#Check for functions
AC_CHECK_FUNCS([fesetround pow rint sqrt])

AC_SUBST(CONF_LDADD, $CONF_LDADD)


#Generate Doxygen file
AM_CONDITIONAL([HAVE_DOXYGEN],
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/Doxyfile])])


# Output files
AC_CONFIG_FILES([Makefile include/Makefile src/Makefile main/Makefile main/analyze/Makefile main/analyze/gparity/Makefile main/tools/Makefile main/fit/Makefile main/fit/gparity/Makefile main/misc/Makefile example/Makefile test/Makefile docs/Makefile])
AC_OUTPUT
AM_INIT_AUTOMAKE([subdir-objects])

