noinst_LIBRARIES = libPathogenesis.a

libPathogenesis_a_SOURCES = PathogenesisModel.cpp Mueller.cpp Pyrogen.cpp Predet.cpp

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 
