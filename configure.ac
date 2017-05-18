#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([include/distribution.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

CPPFLAGS="-Wno-deprecated-declarations $CPPFLAGS"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_REQUIRE_CPP


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
AC_CHECK_LIB([boost_python-py27], [main], [LIBS="-lboost_python-py27 $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_python-py27])])
fi

#Check for boost serialization
AC_CHECK_LIB([boost_serialization], [main], [LIBS="-lboost_serialization $LIBS"],
             [AC_MSG_FAILURE([Could not link to boost_serialization])])



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

# Output files
AC_OUTPUT(Makefile src/Makefile main/Makefile)
AM_INIT_AUTOMAKE([subdir-objects])


