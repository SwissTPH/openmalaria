noinst_PROGRAMS = openMalaria
openMalaria_SOURCES =	global.cpp drug.cpp proteome.cpp infection.cpp withinHostModel.cpp oldWithinHostModel.cpp intervention.cpp \
			summary.cpp event.cpp human.cpp transmissionModel.cpp noVectorControl.cpp vectorControl.cpp population.cpp simulation.cpp setDemo.cpp caseManagement.cpp EntoModel.cpp GSLWrapper.cpp openMalaria.cpp

openMalaria_DEPENDENCIES =  $(top_builddir)/parser/libparser.a

openMalaria_LDADD = $(top_builddir)/parser/libparser.a @LIB_XERC@ @LIB_GSL@ @LIB_GSLCBLAS@ @LIB_BOINC@ @LIB_BOINC_API@ 

AM_CXXFLAGS = @CXXFLAGS@ -static-libgcc -Wall @BOINC_CFLAGS@ -I$(top_srcdir)/include #-Werror

if DARWIN
export MACOSX_DEPLOYMENT_TARGET=10.4
openMalaria_LDFLAGS = @BOINC_LDFLAGS@ -framework Carbon
else
openMalaria_LDFLAGS = @BOINC_LDFLAGS@
endif
