ACLOCAL_AMFLAGS = -I m4
SUBDIRS = xsdcxx model #graphics 

all: xsdcxx model

test: all test-scenarios test-checkpointing

test-scenarios: all
	./test.sh ;

test-checkpointing: all
	./testCheckpointing.sh ;

doxygen: config.doxygen
	-doxygen config.doxygen
