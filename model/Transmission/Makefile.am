# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

noinst_LIBRARIES = libTransmission.a

libTransmission_a_SOURCES = TransmissionModel.cpp NonVector.cpp  Vector.cpp  VectorEmergence.cpp PerHost.cpp VectorSpecies.cpp

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 
