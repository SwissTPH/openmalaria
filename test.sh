#!/bin/sh
# With no arguments, run all scenario*.xml files.
# With arguments A,...,Z, only run scenarioA.xml, ..., scenarioZ.xml

numberOfdllsNeeded=`ldd model/openMalaria | wc -l`
numberOfdllsExpected=4
if [ $numberOfdllsNeeded -ne $numberOfdllsExpected ]
	then echo "you probably don't have all project specific libraries statically linked, please check before BOINC deployment!"
fi

rm test/sandbox/* 2>/dev/null
cp model/openMalaria test/sandbox
strip test/sandbox/openMalaria
cd test/sandbox && cp ../original/* . 2>/dev/null

runScenario() {
  mv scenario$number.xml scenario.xml
  echo "Running scenario$number.xml"
  ./openMalaria  > /dev/null 
  wait
  date
  mv output.txt output$number.txt
  test -e output$number.txt && ../original/compareOutputsFloat.py original$number.txt output$number.txt 1
}

if [ "$#" -eq "0" ]
then
    for file in scenario*.xml
    do
	name=${file%.xml}
  	number=${name#scenario} 
  	runScenario
    done
else
    while [ "$#" -ge "1" ]
    do
    	number=$1
    	runScenario
    	shift
    done
fi
