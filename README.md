
[![Stories in Ready](https://badge.waffle.io/SwissTPH/openmalaria.png?label=ready&title=Ready)](https://waffle.io/SwissTPH/openmalaria)
Open Malaria
============

This git-repository contains source-code for Open Malaria, a simulator program
for studying malaria epidemiology and the impacts of interventions against
malaria.

------

For further documentation, take a look at our
[wiki](https://github.com/SwissTPH/openmalaria/wiki).

Also find our code on zenodo, with DOI-references:
[![DOI](https://zenodo.org/badge/15670/SwissTPH/openmalaria.svg)](https://zenodo.org/badge/latestdoi/15670/SwissTPH/openmalaria)

The stable version of OpenMalaria supports schema
[schema version 34](https://github.com/SwissTPH/openmalaria/wiki/GeneratedSchema34Doc)
and is currently maintained in the
__['master' branch](https://github.com/SwissTPH/openmalaria/tree/master)__.

You can download the latest build here:
[schema-34.0](https://github.com/SwissTPH/openmalaria/releases/tag/schema-34.0)

Status of __[master](https://github.com/SwissTPH/openmalaria/tree/master)__
builds:

[![Linux Build Status](https://travis-ci.org/SwissTPH/openmalaria.svg?branch=master)](https://travis-ci.org/SwissTPH/openmalaria)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/8el77m2gg4aqqnqg/branch/master?svg=true)](https://ci.appveyor.com/project/tph-thuering/openmalaria/branch/master)

Status of __[develop](https://github.com/SwissTPH/openmalaria/tree/develop)__
builds:


[![Build Status](https://travis-ci.org/SwissTPH/openmalaria.svg?branch=develop)](https://travis-ci.org/SwissTPH/openmalaria)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/8el77m2gg4aqqnqg/branch/develop?svg=true)](https://ci.appveyor.com/project/tph-thuering/openmalaria/branch/develop)

Code Coverage:
![Coverage Status](https://coveralls.io/repos/SwissTPH/openmalaria/badge.svg)

License: [GPL v2](http://opensource.org/licenses/GPL-2.0) (see COPYING).

Version
======

The schema version is specified in the following places (all need to be updated
when releasing a new version):

*   DocumentLoader.h
*   schema/scenario.xsd, demography.xsd, etc. (all XSD files without a version number)
*   copy build/schema/scenario_current.xsd to schema/scenario_XX.xsd
*   test/*.xml — update http://openmalaria.org/schema/scenario_XX and (optionally) schemaVersion="XX"
*   version.txt — needed for build service

In theory the "schema namespace version" doesn't need to match the "OpenMalaria"
version and we could update the latter without requiring changes to XML files,
however currently we keep both synchronised (in some ways this is simpler).


Installation instructions:
==================
See [INSTALLATION](https://github.com/SwissTPH/openmalaria/wiki/Installation).

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
