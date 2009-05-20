noinst_LIBRARIES = libxsdcxx.a
libxsdcxx_a_SOURCES = inputData.cpp scenario.cxx
BUILT_SOURCES = scenario.hxx scenario.ixx scenario.cxx

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 

#configure decides whether to use xsd or xsdcxx (must be in PATH)

# use UpperCamelCase type names, java function naming (getFoo(), etc),
#put everything in the namespace scnXml,
#inline simple functions,
#generate doxygen doc,
#VS only: --generate-intellisense
XSDFLAGS = --type-naming ucc --function-naming java\
--namespace-map =scnXml\
--generate-inline\
--generate-doxygen


%.cxx %.hxx %.ixx: %.xsd
	$(XSD) cxx-tree $(XSDFLAGS) $<

clean: clean-am clean-scenario

clean-scenario:
	rm -f scenario.cxx scenario.hxx scenario.ixx
