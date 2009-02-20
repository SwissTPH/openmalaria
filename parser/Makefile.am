noinst_LIBRARIES = libparser.a

libparser_a_SOURCES =	EIRDaily.cpp  documentParameter.cpp  demography.cpp  scenario.cpp  \
			proteome.cpp distribution.cpp content.cpp protein.cpp mutation.cpp proteomeInstance.cpp allele.cpp \
			genoType.cpp ITNDescription.cpp IPTDescription.cpp mass.cpp  ageGroup.cpp  surveys.cpp \
			IRS.cpp  intNode.cpp  inputData.cpp  group.cpp  item.cpp  doubleNode.cpp  \
			vaccineDescription.cpp  ageSpecific.cpp  DOMParser.cpp  intervention.cpp  parameters.cpp \
			interventions.cpp  byAgeItems.cpp  parameter.cpp  CFRGroup.cpp value.cpp  continuous.cpp  \
			entoData.cpp  DOMTreeErrorReporter.cpp  monitoring.cpp  document.cpp  CFR.cpp  timed.cpp  \
			treatmentDetails.cpp  converter.cpp healthSystem.cpp medicate.cpp endPoint.cpp \
			entryPoint.cpp decisions.cpp decision.cpp stringNode.cpp \
			drugRegimen.cpp caseManagement.cpp caseManagements.cpp

AM_CXXFLAGS = @CXXFLAGS@ -I$(top_srcdir)/model -I$(top_srcdir)/include -static-libgcc -Wall @BOINC_CFLAGS@ #-Werror 



