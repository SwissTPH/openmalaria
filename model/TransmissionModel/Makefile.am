noinst_LIBRARIES = libTransmissionModel.a

libTransmissionModel_a_SOURCES = NonVector.cpp  Vector.cpp  VectorInternal.cpp PerHost.cpp VectorSpecies.cpp

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 
