# This is only intended to be a _very simple_ build script for CaseMangementTree.cpp

OM_INCLUDE:=-I../../include
LIBXML_PP_FLAGS:=`pkg-config libxml++-2.6 --cflags --libs`

all : CaseManagementTree

CaseManagementTree : CaseManagementTree.cpp
	c++ -g ${OM_INCLUDE} ${LIBXML_PP_FLAGS} -o CaseManagementTree CaseManagementTree.cpp

clean:
	rm -f CaseManagementTree
