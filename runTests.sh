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

RUN="test"
CMD_PREFIX=""
CMD_MAIN="./openMalaria --scenario"
CMD_POSTFIX=""

runCMD() {
  echo "\033[1;32m"$CMD"\033[0;00m"
  if [ "x$RUN" != "x" ]
  then
    $CMD
  fi
}

runScenario() {
  # delete old checkpoints; necessary after a previous run:
  rm -f checkpoint* seed*
  CMD="$CMD_PREFIX$CMD_MAIN scenario$number.xml $CMD_POSTFIX"
  echo "\033[0;33m"`date`
  touch timeFile
  # First run:
  runCMD
  # While no output, a checkpoint has been written, and $CMD exits successfully:
  while [ ! -f output.txt -a checkpoint -nt timeFile ] && runCMD
  do
    wait
  done
  if [ -f stderr.txt ]
  then
    mv stderr.txt stderr$number.txt
  fi
  if [ -f output.txt ]
  then
    mv output.txt output$number.txt
    echo -n "\033[1;34m"
    ../original/compareOutputsFloat.py original$number.txt output$number.txt 1
  elif [ "$RUN" = "test" ]
  then
    echo "\033[0;31mNo results output; error messages:"
    test -f stderr$number.txt && cat stderr$number.txt
  fi
  echo "\033[0;00m"
}

while [ "$#" -ge "1" ]
do
  if [ "$1" = "--gdb" ]
  then
    CMD_PREFIX="gdb --args "
  elif [ "$1" = "--valgrind" ]
  then
    CMD_PREFIX="valgrind --gen-suppressions=yes --leak-check=full "
  elif [ "$1" = "--valgrind=track" ]
  then
    CMD_PREFIX="valgrind --gen-suppressions=yes --leak-check=full --track-origins=yes "
  elif [ "$1" = "--dont-run" ]
  then
    RUN=""
  elif [ "$1" = "--valid" ]
  then
    CMD_MAIN="xmllint --noout --schema scenario_6.xsd"
    RUN="run"
  elif [ "$1" = "--help" ]
  then
    echo "Usage: \033[1;32m$0 [options] [scenarios]\033[0;00m"
    echo
    echo "Scenarios to be run must be of the form scenarioXX.xml; if any are passed on"
    echo "the command line, XX is substituted for each given; if not then all files of"
    echo "the form scenario*.xml are run as test scenarios."
    echo
    echo "Options:\033[0;33m"
    echo "  --gdb		Run openMalaria through gdb."
    echo "  --valgrind		Run openMalaria through valgrind."
    echo "  --valgrind=track	As --valgrind, but pass --track-origins=yes option (½ performance)."
    echo "  --valid		Validate the XML file(s) using xmllint and the latest schema."
    echo "  --dont-run		Don't actually run openMalaria, just output the commandline."
    echo "  --help		Print this message."
    echo "\033[0;00mOther options starting '--' will be passed to openMalaria. openMalaria options:\033[0;33m"
    ./openMalaria --help
    echo -n "\033[0;00m"
    exit 1;
  elif [ `expr "$1" : '--'` -eq 2 ]
  then
    CMD_POSTFIX="$CMD_POSTFIX $1"
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
