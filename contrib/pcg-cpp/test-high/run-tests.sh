#!/bin/sh
# 
# PCG Random Number Generation for C.
# 
# Copyright 2014-2017 Melissa O'Neill <oneill@pcg-random.org>,
#                     and the PCG Project contributors.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#
# Licensed under the Apache License, Version 2.0 (provided in
# LICENSE-APACHE.txt and at http://www.apache.org/licenses/LICENSE-2.0)
# or under the MIT license (provided in LICENSE-MIT.txt and at
# http://opensource.org/licenses/MIT), at your option. This file may not
# be copied, modified, or distributed except according to those terms.
#
# Distributed on an "AS IS" BASIS, WITHOUT WARRANTY OF ANY KIND, either
# express or implied.  See your chosen license for details.
#
# For additional information about the PCG random number generation scheme,
# visit http://www.pcg-random.org/.
#

echo Performing a quick sanity check...

mkdir -p actual
rm -f actual/*

./check-pcg32 > actual/check-pcg32.out
./check-pcg64 > actual/check-pcg64.out
./check-pcg64_fast > actual/check-pcg64_fast.out

find actual -type f -size -80c -delete

if diff -ru -x .gitignore expected actual
then
    echo All tests succeeded.
else
    echo ''
    if diff -x "*-pcg64_[ck]*.out" \
            -x "*-pcg128_[ck]*.out" -ru expected actual > /dev/null
    then
        echo All tests except tests awkward tests with 128-bit math succceed.
    else
        echo ERROR: Some tests failed.
    fi
fi
