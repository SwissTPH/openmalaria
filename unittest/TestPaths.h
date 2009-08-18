#ifndef Hmod_TestPaths
#define Hmod_TestPaths
// If this file is included unconfigured, this will fail.
// When configured, @CMAKE_CONFIGURED@ is replaced by 1.
#if @CMAKE_CONFIGURED@ != 1
assert (false);
#endif

const char* UnittestSourceDir = "@CMAKE_CURRENT_SOURCE_DIR@/";	// must end with '/'
const char* UnittestScenario = "@CMAKE_CURRENT_SOURCE_DIR@/scenario.xml";
#endif
