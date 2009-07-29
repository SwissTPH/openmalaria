/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
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
CxxTest::StaticSuiteDescription suiteDescription_VectorEmergenceSuite( "VectorEmergenceSuite.h", 50, "VectorEmergenceSuite", suite_VectorEmergenceSuite, Tests_VectorEmergenceSuite );

static class TestDescription_VectorEmergenceSuite_testCalcUpsilonOneHost : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcUpsilonOneHost() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 123, "testCalcUpsilonOneHost" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcUpsilonOneHost(); }
} testDescription_VectorEmergenceSuite_testCalcUpsilonOneHost;

static class TestDescription_VectorEmergenceSuite_testCalcSvDiff : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcSvDiff() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 132, "testCalcSvDiff" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcSvDiff(); }
} testDescription_VectorEmergenceSuite_testCalcSvDiff;

static class TestDescription_VectorEmergenceSuite_testCalcLambda : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcLambda() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 143, "testCalcLambda" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcLambda(); }
} testDescription_VectorEmergenceSuite_testCalcLambda;

static class TestDescription_VectorEmergenceSuite_testCalcXP : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcXP() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 149, "testCalcXP" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcXP(); }
} testDescription_VectorEmergenceSuite_testCalcXP;

static class TestDescription_VectorEmergenceSuite_testCalcPSTS : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcPSTS() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 163, "testCalcPSTS" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcPSTS(); }
} testDescription_VectorEmergenceSuite_testCalcPSTS;

static class TestDescription_VectorEmergenceSuite_testFuncX : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testFuncX() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 173, "testFuncX" ) {}
 void runTest() { suite_VectorEmergenceSuite.testFuncX(); }
} testDescription_VectorEmergenceSuite_testFuncX;

static class TestDescription_VectorEmergenceSuite_testCalcSpectralRadius : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcSpectralRadius() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 185, "testCalcSpectralRadius" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcSpectralRadius(); }
} testDescription_VectorEmergenceSuite_testCalcSpectralRadius;

static class TestDescription_VectorEmergenceSuite_testCalcInv1minusA : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalcInv1minusA() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 189, "testCalcInv1minusA" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalcInv1minusA(); }
} testDescription_VectorEmergenceSuite_testCalcInv1minusA;

static class TestDescription_VectorEmergenceSuite_testCalSvfromEIRdata : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testCalSvfromEIRdata() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 197, "testCalSvfromEIRdata" ) {}
 void runTest() { suite_VectorEmergenceSuite.testCalSvfromEIRdata(); }
} testDescription_VectorEmergenceSuite_testCalSvfromEIRdata;

static class TestDescription_VectorEmergenceSuite_testWholeCalculation : public CxxTest::RealTestDescription {
public:
 TestDescription_VectorEmergenceSuite_testWholeCalculation() : CxxTest::RealTestDescription( Tests_VectorEmergenceSuite, suiteDescription_VectorEmergenceSuite, 204, "testWholeCalculation" ) {}
 void runTest() { suite_VectorEmergenceSuite.testWholeCalculation(); }
} testDescription_VectorEmergenceSuite_testWholeCalculation;

#include <cxxtest/Root.cpp>
