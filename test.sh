#!/bin/sh
# With no arguments, run all scenario*.xml files.
# With arguments A,...,Z, only run scenarioA.xml, ..., scenarioZ.xml

numberOfdllsNeeded=`ldd model/openMalaria | wc -l`
numberOfdllsExpected=4
if [ $numberOfdllsNeeded -ne $numberOfdllsExpected ]
	then echo "you probably don't have all project specific libraries statically linked, please check before BOINC deployment!"
fi

rm -rf test/sandbox 2>/dev/null
mkdir -p test/sandbox
cp model/openMalaria test/sandbox/openMalaria
strip test/sandbox/openMalaria
cd test/sandbox && cp ../original/* . 2>/dev/null

runScenario() {
  mv scenario$number.xml scenario.xml
  echo "Running scenario$number.xml"
  ./openMalaria  > /dev/null 
  wait
  date
  if [ -e stderr.txt ]
  then
    mv stderr.txt stderr$number.txt
  fi
  if [ -e output.txt ]
    then
      mv output.txt output$number.txt
      ../original/compareOutputsFloat.py original$number.txt output$number.txt 1
    else
      echo "No output; error messages:"
      cat stderr$number.txt
  fi
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
