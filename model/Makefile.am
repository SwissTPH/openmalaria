# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

SUBDIRS = TransmissionModel WithinHostModel Drug Pathogenesis Clinical
noinst_PROGRAMS = openMalaria
openMalaria_SOURCES = openMalaria.cpp\
 GSLWrapper.cpp timer.cpp BoincWrapper.cpp\
 global.cpp summary.cpp population.cpp simulation.cpp \
 proteome.cpp\
 WithinHostModel.cpp\
 intervention.cpp human.cpp\
 TransmissionModel.cpp EntoIntervention.cpp\
 InfectionIncidenceModel.cpp\
 gzstream.C

openMalaria_DEPENDENCIES =  $(top_builddir)/xsdcxx/libxsdcxx.a\
	TransmissionModel/libTransmissionModel.a\
	WithinHostModel/libWithinHostModel.a Drug/libDrug.a\
	Pathogenesis/libPathogenesis.a\
	Clinical/libClinical.a

openMalaria_LDADD = $(openMalaria_DEPENDENCIES)\
 @LIB_PATH@ @LIB_BOINC@ @LIB_BOINC_API@ @LIB_STATIC@

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I$(top_srcdir)/include -I@top_builddir@/xsdcxx
AM_CXXFLAGS = @CXXFLAGS@ -Wall #-Werror

if DARWIN
export MACOSX_DEPLOYMENT_TARGET=10.4
openMalaria_LDFLAGS = @BOINC_LDFLAGS@ -framework Carbon
else
openMalaria_LDFLAGS = @BOINC_LDFLAGS@
endif
