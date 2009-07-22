#!/bin/sh
# With no arguments, run all scenario*.xml files.
# With arguments A,...,Z, only run scenarioA.xml, ..., scenarioZ.xml

# replaced by CMake
TEST_SOURCE=@CMAKE_CURRENT_SOURCE_DIR@
TEST_BINARY=@CMAKE_CURRENT_BINARY_DIR@
cd $TEST_BINARY
if [ ! -d $TEST_SOURCE ]
then
  echo "Don't run this script directly; configure CMake then use the version in the CMake build dir."
  exit 1;
fi
# if same dir:
if [ $TEST_SOURCE -ef $TEST_BINARY ]
then
  echo "This must not be run from the source dir!"
  exit 1
fi

# executable
OM_NAME="openMalaria"
OM_BIN="../$OM_NAME"
# cmake 2.4 still puts it in the model dir, I think:
if [ ! -x $OM_BIN -o ../model/$OM_NAME -nt $OM_BIN ]
then
  OM_BIN="../model/$OM_NAME"
fi
if [ ! -x $OM_BIN ]
then
  echo "Not found: $OM_NAME. Please compile."
  exit 1
fi

# If RUN="test", then an error is printed when the output.txt file doesn't exist:
RUN="test"
# Used to run through gdb/valgrind:
CMD_PREFIX=""
# The command to run (plus any args before the scenario file):
CMD_MAIN="$OM_BIN --scenario"
# Extra arguments to pass to the command:
CMD_POSTFIX=""
# If true, don't print so many messages:
QUIET=""

# Echo the command line to be run and run it
runCMD() {
  if [ "x$QUIET" = "x" ]
  then
    echo "\033[1;32m"$CMD"\033[0;00m"
  fi
  if [ "x$RUN" != "x" ]
  then
    $CMD
  fi
}

# Run, with file scenario$number.xml
runScenario() {
  # delete old checkpoints; necessary after a previous run:
  # Also delete the output.txt file if existing since it won't get overridden now
  rm -f checkpoint* seed* output.txt
  CMD="$CMD_PREFIX$CMD_MAIN $TEST_SOURCE/scenario$number.xml $CMD_POSTFIX"
  if [ "x$QUIET" = "x" ]
  then
    echo "\033[0;33m"`date`
  fi
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
    $TEST_SOURCE/compareOutputsFloat.py $TEST_SOURCE/original$number.txt output$number.txt 1
  elif [ "$RUN" = "test" ]
  then
    echo "\033[0;31mNo results output; error messages:"
    test -f stderr$number.txt && cat stderr$number.txt
  fi
  echo "\033[0;00m"
}

# Test for options
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
  elif [ "$1" = "--quiet" ]
  then
    QUIET="1>/dev/null 2>/dev/null"
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

# Run all tests:
if [ "$#" -eq "0" ]
then
  if [ ! -f $TEST_SOURCE/scenario2.xml ] # make sure at least one exists
    then
      echo "Error: no scenario files present (at least not scenario2.xml)"
      exit 1
    fi
    for path in $TEST_SOURCE/scenario*.xml
    do
	file=`basename $path`
	name=${file%.xml}
  	number=${name#scenario} 
  	runScenario
    done
# or selected tests:
else
    while [ "$#" -ge "1" ]
    do
    	number=$1
    	runScenario
    	shift
    done
fi
