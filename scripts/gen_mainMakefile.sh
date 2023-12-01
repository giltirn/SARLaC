#!/usr/bin/bash

INST_DIR=$1

#From a directory containing .C main programs, generate a Makefile.am

exec> Makefile.am

echo "exampledir = \$(prefix)/example"
echo "testdir = \$(prefix)/test"

echo -n $INST_DIR"_PROGRAMS ="
for i in $(ls *.C); do
    STUB=$(echo $i | sed 's/\.C$//')
    echo -n " "$STUB
done
echo ""

for i in $(ls *.C); do
    STUB=$(echo $i | sed 's/\.C$//')
    echo $STUB"_SOURCES = "$i
done

echo "AM_CPPFLAGS = -I\$(top_srcdir)/include -I\$(srcdir)/include"
echo "LDADD =\$(top_builddir)/src/libcpsfit.a"
if [ -e subdirs.inc ]; then cat subdirs.inc ; fi
