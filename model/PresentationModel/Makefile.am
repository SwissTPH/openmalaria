noinst_LIBRARIES = libPresentationModel.a

libPresentationModel_a_SOURCES =  Mueller.cpp Pyrogen.cpp Predet.cpp

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 

clean: clean-am clean-scenario

clean-scenario:
	rm -f scenario.cxx scenario.hxx scenario.ixx
