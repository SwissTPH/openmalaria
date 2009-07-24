/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ParenPrinter.h>

int main() {
 return CxxTest::ParenPrinter().run();
}
#include "PerHostSuite.h"

static PerHostSuite suite_PerHostSuite;

static CxxTest::List Tests_PerHostSuite = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_PerHostSuite( "PerHostSuite.h", 28, "PerHostSuite", suite_PerHostSuite, Tests_PerHostSuite );

static class TestDescription_PerHostSuite_testRelativeAvailability : public CxxTest::RealTestDescription {
public:
 TestDescription_PerHostSuite_testRelativeAvailability() : CxxTest::RealTestDescription( Tests_PerHostSuite, suiteDescription_PerHostSuite, 35, "testRelativeAvailability" ) {}
 void runTest() { suite_PerHostSuite.testRelativeAvailability(); }
} testDescription_PerHostSuite_testRelativeAvailability;

#include "VectorEmergenceSuite.h"

static VectorEmergenceSuite suite_VectorEmergenceSuite;

static CxxTest::List Tests_VectorEmergenceSuite = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_VectorEmergenceSuite( "VectorEmergenceSuite.h", 31, "VectorEmergenceSuite", suite_VectorEmergenceSuite, Tests_VectorEmergenceSuite );

static class TestDescription_VectorEmergenceSuite_testDummy : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testDummy() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 40, "testDummy" ) {}
 void runTest() { suite_VectorEmergenceSuite.testDummy(); }
} testDescription_VectorEmergenceSuite_testDummy;

static class TestDescription_VectorEmergenceSuite_testWholeCalculation : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testWholeCalculation() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 44, "testWholeCalculation" ) {}
 void runTest() { suite_VectorEmergenceSuite.testWholeCalculation(); }
} testDescription_VectorEmergenceSuite_testWholeCalculation;

#include <cxxtest/Root.cpp>
