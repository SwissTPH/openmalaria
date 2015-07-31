#!/bin/bash

# This script is used in the build process with travis-ci
#It creates checksums (md5 and sha256) from openMalaria, outputs it to the console and creates files containing the checksum
#It takes the build path of the openMalaria binary and an output directory as positional parameters. If none provided takes the current dir.

BINARY_DIR=$1
OUTPUT_DIR=$2
FILE=openMalaria
echo $PWD

[ -z "$BINARY_DIR" ] && BINARY_DIR=.
[ -z "$OUTPUT_DIR" ] && OUTPUT_DIR=.

if [ "$TRAVIS_OS_NAME" == "linux" ]
then
  md5sum $BINARY_DIR/$FILE | tee $OUTPUT_DIR/$FILE.md5
  sha256sum $BINARY_DIR/$FILE | tee $OUTPUT_DIR/$FILE.sha256
elif [ "$TRAVIS_OS_NAME" == "osx" ]
then
  md5 $BINARY_DIR/$FILE | tee $OUTPUT_DIR/$FILE.md5
  /usr/local/bin/gsha256sum $BINARY_DIR/$FILE | tee $OUTPUT_DIR/$FILE.sha256
fi

exit $?
