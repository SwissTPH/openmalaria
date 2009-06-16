# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

noinst_LIBRARIES = libClinical.a

libClinical_a_SOURCES = ClinicalModel.cpp\
	ImmediateOutcomes.cpp OldCaseManagement.cpp event.cpp\
	EventScheduler.cpp

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I@top_builddir@/xsdcxx -I$(top_srcdir)/include
AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall #-Werror 
