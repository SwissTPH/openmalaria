Open Malaria
============

This git-repository contains source-code for Open Malaria, a simulator program for studying malaria epidemiology and the impacts of interventions against malaria.
------

For further documentation, take a look at our [wiki](https://github.com/SwissTPH/openmalaria/wiki).

The stable version of OpenMalaria supports schema [schema version 33](https://github.com/SwissTPH/openmalaria/wiki/GeneratedSchema33Doc) and is currently maintained in the __['master' branch](https://github.com/SwissTPH/openmalaria/tree/master)__.


Status of __[master](https://github.com/SwissTPH/openmalaria/tree/master)__ builds:

[![Linux Build Status](https://travis-ci.org/SwissTPH/openmalaria.svg?branch=master)](https://travis-ci.org/SwissTPH/openmalaria) [![Windows Build Status](https://ci.appveyor.com/api/projects/status/8el77m2gg4aqqnqg/branch/master?svg=true)](https://ci.appveyor.com/project/tph-thuering/openmalaria/branch/master)

The development version of OpenMalaria supports the upcoming [schema version 34](https://github.com/SwissTPH/openmalaria/wiki/GeneratedSchema34Doc) and is currently maintained in the __['develop' branch](https://github.com/SwissTPH/openmalaria/tree/develop)__.

Status of __[develop](https://github.com/SwissTPH/openmalaria/tree/develop)__ builds:


[![Build Status](https://travis-ci.org/SwissTPH/openmalaria.svg?branch=develop)](https://travis-ci.org/SwissTPH/openmalaria)
![Coverage Status](https://coveralls.io/repos/SwissTPH/openmalaria/badge.svg)


License: [GPL v2](http://opensource.org/licenses/GPL-2.0) (see COPYING).


Build instructions:
===================

```
mkdir build && cd build
ccmake ..
Press 'c', look over options, press 'c' again and 'g'
make -j4
ctest -j4
```

For testing and development, ideally use debug builds (which enable some
asserts to do with simulation time usage).

Code subdirectories:
=============
|- dir    -|- description -|
|----------|:-------------------------------------------------------------------------:|
| contrib | Third-party libraries, distributed under the same repo for convenience.   |
| model   | Source code for the malaria model.                                        |
| test    | High-level testing: test scenarios with expected outputs. Also run-time files: densities.csv, scenario_?.xsd, Nv0scenario*.txt. |
| unittest| Low-level testing: unittests for the model using cxxunit. |
| util    | Extra scripts associated with OpenMalaria. |
| schema  | scenario schema files (see schema/policy.txt for details) |
| schema/scenario.xsd | The latest (partial) schema file. |
| schema/entomology.xsd, schema/demography.xsd, etc | components of the latest schema, included from scenario.xsd. |
| schema/scenario_*.xsd | Copies of released schema versions, with all components inlined in the same file. |

This git repository is currently maintained by members of the [Dynamical Modelling Group](http://www.swisstph.ch/about-us/departments/epidemiology-and-public-health-eph/health-systems-research-and-dynamical-modelling/dynamical-modelling.html) of the __Swiss Tropical and Public Health institute__ and other collaborators.
