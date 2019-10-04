# 
# PCG Random Number Generation for C++.
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

TARGETS           = check-pcg32 check-pcg64 check-pcg64_fast

# special flags for some compilers
CPPFLAGS_clang++ += \
	-Weverything \
	-Wno-unknown-warning-option \
	-Wno-c++98-compat \
	-Wno-c++98-compat-pedantic \
	-Wno-date-time \
	-Wno-undef \
	-Wno-header-hygiene \
	-Wno-unused-macros

CPPFLAGS += -I../include -Wall -Wextra $(CPPFLAGS_$(CXX))
CXXFLAGS += -std=c++11
CC       = $(CXX)       # Cheat so that linking uses the C++ compiler

all: $(TARGETS)

test: $(TARGETS)
	sh run-tests.sh

clean:
	rm -f *.o $(TARGETS)
	rm -rf actual

pcg-test-noadvance.o: pcg-test-noadvance.cpp ../include/pcg_random.hpp \
  ../include/pcg_extras.hpp ../include/pcg_uint128.hpp
pcg-test.o: pcg-test.cpp ../include/pcg_random.hpp \
  ../include/pcg_extras.hpp ../include/pcg_uint128.hpp
check-pcg32.o: check-pcg32.cpp pcg-test.cpp ../include/pcg_random.hpp \
  ../include/pcg_extras.hpp ../include/pcg_uint128.hpp
check-pcg64.o: check-pcg64.cpp pcg-test.cpp ../include/pcg_random.hpp \
  ../include/pcg_extras.hpp ../include/pcg_uint128.hpp
check-pcg64_fast.o: check-pcg64_fast.cpp pcg-test.cpp \
  ../include/pcg_random.hpp ../include/pcg_extras.hpp \
  ../include/pcg_uint128.hpp
