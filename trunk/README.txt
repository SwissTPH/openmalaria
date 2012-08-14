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


Scenario schema files (see schema/policy.txt for details):
schema/scenario.xsd : The latest (partial) schema file.
schema/entomology.xsd, schema/demography.xsd, etc: components of the latest schema, included from scenario.xsd.
schema/scenario_*.xsd : Copies of released schema versions, with all components inlined in the same file.
