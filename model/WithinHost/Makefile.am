# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

noinst_LIBRARIES = libWithinHost.a

libWithinHost_a_SOURCES = WithinHostModel.cpp\
  Descriptive.cpp OldIPT.cpp Dummy.cpp\
  Infection.cpp DescriptiveInfection.cpp DummyInfection.cpp OldIPTInfection.cpp\
  Empirical.cpp EmpiricalInfection.cpp

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 
