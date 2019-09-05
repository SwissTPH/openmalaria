This directory contains XSD schema files describing scenario files.

The latest schema is composed of the unversion files (scenario.xsd,
demography.xsd, etc.). The schema is compiled into a single file every time
OpenMalaria is built: build/schema/scenario_current.xsd, since a single file is
easier to use when handling validations.

Update policy:
*   update scenario.xsd and C++ code
*   update the schema translator [translateXML.py](https://github.com/vecnet/openmalaria.tools/blob/master/openmalaria/tools/translateXML.py) with a function translating the last release version to the next (translate_X_to_Y where X is current version, Y=X+1), doing necessary translations and writing a comment explaining the change
*   update test scenarios when necessary (at least those run by ctest)
*   in DocumentLoader.h, update SCHEMA_VERSION to the next version number (otherwise it'll refuse to run updated scenarios)
*   update `version.txt`

Release policy (new schema versions):
*   copy build/schema/scenario_current.xsd to scenario_VER.xsd
*   update the CURRENT_VERSION var in the schema translator
*   run all scenarios in ../test dir through the schema translator (those already updated by hand should have had their version numbers updated already, so the schema translator doesn't try updating them again)
