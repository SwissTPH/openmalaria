#!/bin/bash

if [ -z "$TRAVIS_BUILD_ID" ]
then
	echo "Error: Use only for automated Travis-ci builds"  1>&2
	exit 1
fi

BUILD=$BASE_URL'/builds/'$TRAVIS_BUILD_ID
COMMIT=$BASE_URL'/commits/'$TRAVIS_BRANCH'/'$TRAVIS_COMMIT
JOB=$BASE_URL'/jobs/'$TRAVIS_JOB_ID
OS=$TRAVIS_OS_NAME
TAG=$TRAVIS_TAG

echo -e '{ "'$TRAVIS_JOB_ID'" : {' \
      '"build" : "'$BUILD'", '"\n" \
      '"commit" : "'$COMMIT'", '"\n" \
      '"job" : "'$JOB'",'"\n" \
      '"os" : "'$OS'",'"\n" \
      '"tag" : "'$TAG'"'"\n" \
      '}}'
