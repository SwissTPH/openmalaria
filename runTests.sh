#!/bin/sh
# With no arguments, run all scenario*.xml files.
# With arguments A,...,Z, only run scenarioA.xml, ..., scenarioZ.xml

# ldd test disabled: 4 isn't always the right number; this information isn't useful to most people running the tests.
#numberOfdllsNeeded=`ldd model/openMalaria | wc -l`
#numberOfdllsExpected=4
#if [ $numberOfdllsNeeded -ne $numberOfdllsExpected ]
#	then echo "you probably don't have all project specific libraries statically linked, please check before BOINC deployment!"
#fi

rm -rf test/sandbox/* 2>/dev/null
mkdir -p test/sandbox
cp model/openMalaria test/sandbox/openMalaria
#strip test/sandbox/openMalaria
cd test/sandbox && cp ../original/* . 2>/dev/null

CHECKPOINT=""
PRINT_MODEL=""
CMD_PREFIX=""

runScenario() {
  # delete old checkpoints; necessary after a previous run:
  rm -f checkpoint* seed*
  CMD="$CMD_PREFIX./openMalaria --scenario scenario$number.xml $CHECKPOINT $PRINT_MODEL"
  date
  echo $CMD
  touch timeFile
  $CMD
  wait
  # Checkpoint written after timeFile was last touched, and no output:
  while [ ! -f output.txt -a checkpoint -nt timeFile ]
  do
    echo $CMD
    $CMD
    wait
  done
  if [ -e stderr.txt ]
  then
    mv stderr.txt stderr$number.txt
  fi
  if [ -e output.txt ]
    then
      mv output.txt output$number.txt
      ../original/compareOutputsFloat.py original$number.txt output$number.txt 1
    else
      echo "No results output; error messages:"
      cat stderr$number.txt
  fi
  echo
}

while [ "$#" -ge "1" ]
do
  if [ "$1" = "--checkpoint" ]
  then
    CHECKPOINT=$1
  elif [ "$1" = "--gdb" ]
  then
    CMD_PREFIX="gdb --args "
  elif [ "$1" = "--print-model" ]
  then
    PRINT_MODEL=$1
  elif [ "$1" = "--help" ]
  then
    echo "$0 [options] [scenario numbers]"
    echo "If no scenario numbers are given, all scenarios are run."
    echo "Options:"
    echo "  --checkpoint	Run checkpointing tests. These no longer require BOINC."
    echo "  --gdb		Run openMalaria through gdb"
    echo "  --help		Print this message."
    echo "Other options supported by openMalaria will be passed."
    exit 1;
  else
    break
  fi
  shift
done


if [ "$#" -eq "0" ]
then
    if [ ! -f scenario2.xml ] # make sure at least one exists
    then
      echo "Error: no scenario files present (at least not scenario2.xml)"
      exit 1
    fi
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
