noinst_PROGRAMS = openMalaria
openMalaria_SOURCES = openMalaria.cpp\
 GSLWrapper.cpp timer.cpp boincWrapper.cpp\
 global.cpp summary.cpp population.cpp simulation.cpp \
 drug.cpp proteome.cpp\
 Infection.cpp DescriptiveInfection.cpp DummyInfection.cpp\
 withinHostModel.cpp oldWithinHostModel.cpp OldIPTWithinHostModel.cpp DummyWithinHostModel.cpp\
 intervention.cpp event.cpp human.cpp\
 transmissionModel.cpp noVectorControl.cpp VectorControlInternal.cpp vectorControl.cpp\
 CaseManagementModel.cpp OldCaseManagement.cpp NewCaseManagement.cpp\
 EntoIntervention.cpp\
 muellerMorbidityModel.cpp morbidityModel.cpp predetMorbidityModel.cpp pyrogenMorbidityModel.cpp

openMalaria_DEPENDENCIES =  $(top_builddir)/xsdcxx/libxsdcxx.a

openMalaria_LDADD = $(top_builddir)/xsdcxx/libxsdcxx.a -L../lib -L$BOINC_LIB_PATH @LIB_XERC@ @LIB_GSL@ @LIB_GSLCBLAS@ @LIB_BOINC@ #@LIB_BOINC_API@ 

AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall @BOINC_CFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/xsdcxx #-Werror

if DARWIN
export MACOSX_DEPLOYMENT_TARGET=10.4
openMalaria_LDFLAGS = @BOINC_LDFLAGS@ -framework Carbon
else
openMalaria_LDFLAGS = @BOINC_LDFLAGS@
endif
