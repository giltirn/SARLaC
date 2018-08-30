if HAVE_DOXYGEN
#Based on instructions from http://chris-miceli.blogspot.com/2011/01/integrating-doxygen-with-autotools.html

doxyfile.stamp:
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp

html-local: doxyfile.stamp
pdf-local: doxyfile.stamp
	$(MAKE) -C latex
	cp -s latex/refman.pdf .

clean-local:
	rm -rf $(top_builddir)/docs/html
	rm -rf $(top_builddir)/docs/latex
endif 
