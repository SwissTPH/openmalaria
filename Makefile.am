# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

ACLOCAL_AMFLAGS = -I m4
SUBDIRS = xsdcxx model #graphics 

all: xsdcxx model
	mkdir -p test && test -d test/original || ln -s @abs_top_srcdir@/test/original test/
	test -x runTests.sh || ln -s @top_srcdir@/runTests.sh ./

test: all
	./runTests.sh --checkpoint

doxygen: config.doxygen
	-doxygen config.doxygen
