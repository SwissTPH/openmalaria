noinst_PROGRAMS = openMalaria
openMalaria_SOURCES =	global.cpp drug.cpp proteome.cpp Infection.cpp DescriptiveInfection.cpp withinHostModel.cpp oldWithinHostModel.cpp intervention.cpp \
			summary.cpp event.cpp human.cpp transmissionModel.cpp noVectorControl.cpp VectorControlInternal.cpp vectorControl.cpp population.cpp simulation.cpp caseManagement.cpp EntoIntervention.cpp GSLWrapper.cpp openMalaria.cpp timer.cpp

openMalaria_DEPENDENCIES =  $(top_builddir)/xsdcxx/libxsdcxx.a

openMalaria_LDADD = $(top_builddir)/xsdcxx/libxsdcxx.a @LIB_XERC@ @LIB_GSL@ @LIB_GSLCBLAS@ @LIB_BOINC@ @LIB_BOINC_API@ 

AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall @BOINC_CFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/xsdcxx #-Werror

if DARWIN
export MACOSX_DEPLOYMENT_TARGET=10.4
openMalaria_LDFLAGS = @BOINC_LDFLAGS@ -framework Carbon
else
openMalaria_LDFLAGS = @BOINC_LDFLAGS@
endif
