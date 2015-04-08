Open Malaria
============

[![Build Status](https://travis-ci.org/SwissTPH/openmalaria.svg)](https://travis-ci.org/SwissTPH/openmalaria)

License: [GPL v2](http://opensource.org/licenses/GPL-2.0) (see COPYING).


Build instructions:
===================

`mkdir build && cd build`
`ccmake ..`
`Press 'c', look over options, press 'c' again and 'g'`
`make -j4`
`ctest -j4`

For testing and development, ideally use debug builds (which enable some
asserts to do with simulation time usage).

For further documentation, see the wiki:
[https://github.com/SwissTPH/openmalaria/wiki/Start](https://github.com/SwissTPH/openmalaria/wiki/Start)


Code subdirs:
=============
contrib	Third-party libraries, distributed under the same repo for convenience.
include		Header files associated with model.
model		Source code for the malaria model.
test		High-level testing: test scenarios with expected outputs. Also run-time files: densities.csv, scenario_?.xsd, Nv0scenario*.txt.
unittest	Low-level testing: unittests for the model using cxxunit.
util		Extra scripts associated with OpenMalaria.


Scenario schema files (see schema/policy.txt for details):
schema/scenario.xsd : The latest (partial) schema file.
schema/entomology.xsd, schema/demography.xsd, etc: components of the latest schema, included from scenario.xsd.
schema/scenario_*.xsd : Copies of released schema versions, with all components inlined in the same file.
