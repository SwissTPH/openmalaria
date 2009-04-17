ACLOCAL_AMFLAGS = -I m4
SUBDIRS = xsdcxx model #graphics 

all: xsdcxx model
	mkdir -p test && test -d test/original || ln -s @abs_top_srcdir@/test/original test/
	test -x runTests.sh || ln -s @top_srcdir@/runTests.sh ./

test: all
	./runTests.sh

doxygen: config.doxygen
	-doxygen config.doxygen
