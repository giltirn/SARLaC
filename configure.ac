#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([include/distribution.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([PYTHON], [python-2.7])


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

AC_PROG_RANLIB

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_OUTPUT(Makefile src/Makefile main/Makefile)
AM_INIT_AUTOMAKE([subdir-objects])
