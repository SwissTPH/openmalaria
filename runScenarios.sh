#!/bin/sh
# Runs all scenarios with arguments passed
#doesn't check output

rm -rf test/sandbox/* 2>/dev/null
mkdir -p test/sandbox
cp model/openMalaria test/sandbox/openMalaria
strip test/sandbox/openMalaria
cd test/sandbox && cp ../original/* . 2>/dev/null

    for file in scenario*.xml
    do
	name=${file%.xml}
  	number=${name#scenario} 
  mv scenario$number.xml scenario.xml
  echo "Running scenario$number.xml"
  ./openMalaria "$*"
  wait
    done
