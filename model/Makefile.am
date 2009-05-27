SUBDIRS = Pathogenesis TransmissionModel WithinHostModel
noinst_PROGRAMS = openMalaria
openMalaria_SOURCES = openMalaria.cpp\
 GSLWrapper.cpp timer.cpp boincWrapper.cpp\
 global.cpp summary.cpp population.cpp simulation.cpp \
 drug.cpp proteome.cpp\
 WithinHostModel.cpp\
 intervention.cpp event.cpp human.cpp\
 TransmissionModel.cpp\
 CaseManagementModel.cpp OldCaseManagement.cpp NewCaseManagement.cpp\
 EntoIntervention.cpp\
 gzstream.C

openMalaria_DEPENDENCIES =  $(top_builddir)/xsdcxx/libxsdcxx.a Pathogenesis/libPathogenesis.a WithinHostModel/libWithinHostModel.a TransmissionModel/libTransmissionModel.a

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
