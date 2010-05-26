Open Malaria

Licence: GPL v2 (see COPYING).


Brief build instructions:

mkdir build && cd build
ccmake ..
Press 'c', look over options, press 'c' again and 'g'
make
make test

For further documentation, see the wiki:
http://code.google.com/p/openmalaria/wiki/Start


Code subdirs:
contrib	Third-party libraries, distributed under the same repo for convenience.
graphics	Separate graphics app. Currently out-of-date (untested with latest simulator code).
include		Header files associated with model & xsdcxx code.
model		Source code for the malaria model.
test		High-level testing: test scenarios with expected outputs. Also run-time files: densities.csv, scenario_?.xsd, Nv0scenario*.txt.
unittest	Low-level testing: unittests for the model using cxxunit.
util		Extra scripts associated with OpenMalaria.


Scenario schema files:
model/scenario.xsd	The latest schema file. Code is generated from this file.
test/scenario_*.xsd	Copies of all schema versions. The latest version here should _always_ be a copy of model/scenario.xsd. Scenario XML files refer to one of these schemas, used for validation (usually the latest schema version, but not required to be).
