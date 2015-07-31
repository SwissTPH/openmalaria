#!/bin/bash

echo -e "[$TRAVIS_JOB_ID : {" \
      "build: $BASE_URL/builds/$TRAVIS_BUILD_ID,\n" \
      "commit: $BASE_URL/commits/$TRAVIS_BRANCH/$TRAVIS_COMMIT,\n" \
      "job: $BASE_URL/jobs/$TRAVIS_JOB_ID,\n" \
      "os: $TRAVIS_OS_NAME,\n" \
      "tag: $TRAVIS_TAG\n" \
      "}]"
