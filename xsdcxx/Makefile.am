noinst_LIBRARIES = libxsdcxx.a

libxsdcxx_a_SOURCES = inputData.cpp scenario.cxx

AM_CXXFLAGS = @CXXFLAGS@ -I$(top_srcdir)/xsdcxx -I$(top_srcdir)/include -static-libgcc -Wall @BOINC_CFLAGS@ #-Werror 

#configure decides whether to use xsd or xsdcxx (must be in PATH)

XSDFLAGS = --type-naming ucc --function-naming java --parser-regex '/(.+)/parse\u$$1/' --generate-doxygen

all : inputData.o scenario.o

inputData.o : inputData.cpp scenario.hxx
scenario.o : scenario.cxx scenario.hxx

%.cxx %.hxx %.ixx: ../test/original/%.xsd
	$(XSD) cxx-tree $(XSDFLAGS) $<

clean: clean-am clean-scenario

clean-scenario:
	rm -f scenario.cxx scenario.hxx
