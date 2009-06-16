# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

ACLOCAL_AMFLAGS = -I m4
SUBDIRS = xsdcxx model #graphics 

all: xsdcxx model
	mkdir -p test && test -d test/original || ln -s @abs_top_srcdir@/test/original test/
	test -x runTests.sh || ln -s @top_srcdir@/runTests.sh ./

# Run all valid tests (there may also exist broken tests)
test: all
	./runTests.sh --checkpoint 1 2 3 4 5 6 7 9 11 12 IPT

doxygen: config.doxygen
	-doxygen config.doxygen
