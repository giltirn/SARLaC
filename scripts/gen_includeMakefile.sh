#Generate a Makefile from within the include directory

echo 'SUBDIRS =' > Makefile.am
#Put config.h in main include directory
echo 'include_HEADERS = ../config.h' >> Makefile.am

#Maintain subdirectories
echo -n 'nobase_include_HEADERS = ' >> Makefile.am

for i in $(find . -name '*.h' | sed 's/^\.\///'); do
    echo -n "$i " >> Makefile.am
done
