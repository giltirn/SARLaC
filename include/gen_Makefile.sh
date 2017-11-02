echo 'SUBDIRS =' > Makefile.am
echo 'include_HEADERS = ../config.h' >> Makefile.am
echo -n 'nobase_include_HEADERS = ' >> Makefile.am
for i in $(ls *.h); do
    echo -n "$i " >> Makefile.am
done
