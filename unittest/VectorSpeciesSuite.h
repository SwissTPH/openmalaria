/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef Hmod_VectorSpeciesSuite
#define Hmod_VectorSpeciesSuite

#include "global.h"
#include "configured/TestPaths.h"	// from config; but must be included from the build dir
#include <cxxtest/TestSuite.h>
#include "ExtraAsserts.h"
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <list>
#include <xsd/cxx/tree/exceptions.hxx>

#include "population.h"
#include "human.h"
#include "Transmission/Vector.h"
#include "Transmission/VectorSpecies.h"


/** Read an array into a vector. */
template<class T>
void operator>> (const YAML::Node& node, vector<T>& vec) {
  vec.resize (node.size());
  size_t i = 0;
  for (YAML::Iterator it = node.begin(); it != node.end(); ++it, ++i) {
    vec[i] = it->Read<double>();
  }
}
/** As above, but with function syntax. */
template<class T>
vector<T> yaml2Std (const YAML::Node& node) {
  vector<T> ret;
  node >> ret;
  return ret;
}

/** Read a YAML node into a WeibullDecayedValue.
 *
 * Should be an operator= method of WeibullDecayedValue, but that would make
 * the main simulator depend on yaml-cpp. */
WeibullDecayedValue yaml2WeibullDecayedValue (const YAML::Node& node) {
  double initial, hl, k = 1.0;
  node["initial"] >> initial;
  node["Halflife"] >> hl;
  const YAML::Node *Wbk = node.FindValue ("Weibullk");
  if (Wbk != NULL)
    *Wbk >> k;
  WeibullDecayedValue ret;
  ret.setParameters (initial, hl, k);
  return ret;
}


// The assert we use in calculateEIR unittests. It tests vtm and species params against node.
#   define TS_ASSERT_SPECIES_APPROX(node) VectorSpeciesSuite::doAssertSpecies(__FILE__,__LINE__, (node))


/** Unit test for VectorTransmissionSpecies code.
 *
 * Possible tests to write:
 * 
VectorSpecies unit testing:
  initialise:
    read parameters from XML
    calculate EIR array from fourier series, rotate and add to return param
    read emergence params from XML and scale
  initFeedingCycleProbs:
    calculates P_A, P_df, P_dif
    return total pop avail
  initMainSimulation:
    runs VectorEmergence code
    validates emerge params
    prints out new scenario file
  advancePeriod:
    calculates partialEIR, and some intermediaries
    test with no interventions, and each intervention singly
    depends on: a human population with PerHostTransmission elements
  convertLengthToFullYear:
    resizes an array
  rotateArray:
    offsets all element indicies by some x mod length
 */
class VectorSpeciesSuite : public CxxTest::TestSuite
{
public:
  VectorSpeciesSuite() {
    Global::clResourcePath = UnittestSourceDir;
    ifstream file(Global::lookupResource ("VectorSpeciesSuite.yaml").c_str());
    
    YAML::Parser parser(file);
    parser.GetNextDocument(doc);
    
    file.close();
  }
  
  /** Initialise species.
   *
   * Rather than directly initialise the elements we want to use, try to set up
   * the whole simulation. I thought it'd be easier...
   *
   * Note that using doc["key"] to access map elements isn't speed-optimal when
   * many key-value pairs exist, but for a small number it should be fine. */
  void setUp () {
    try {
    // "unused" node, but it still checks element exists (confirming correct data file)
    /*const YAML::Node& node =*/ doc["VectorSpeciesSuite"];
    
    createDocument (UnittestScenario);
    Global::initGlobal();
    simulation = new Simulation();
    
    // Normally done by Simulation::start():
    simulation->simulationTime = 0;
    simulation->_population->estimateRemovalRates();
    simulation->_population->initialiseHumanList();
    simulation->_population->setupPyramid(false);
    simulation->simulationTime = simulationTime = 1;
    
    // No eir data is present in scenario.xml, so it should directly initialise all model parameters
    vtm = dynamic_cast<VectorTransmission*> (simulation->_population->_transmissionModel);
    assert (vtm != NULL);
    assert (vtm->numSpecies == 1);
    species = &vtm->species[0];
    population = &simulation->_population->_population;
    
    } catch (const ::xsd::cxx::tree::exception<char>& e) {
      cout << "\nXSD Exception: " << e.what() << '\n' << e << endl;
      throw;
    } catch (const exception& e) {
      cout << "\nException: " << e.what() << endl;
      throw;
    } catch (...) {
      cout << "\nUnknown exception" << endl;
      throw;
    }
  }
  
  void tearDown () {
    delete simulation;
  }
  
  void testCalcInverseDFTExp () {
    const YAML::Node& node = doc["calcInverseDFTExp"];
    vector<double> out(10), fc, expected;
    node["fourierCoefficients"] >> fc;
    node["output"] >> expected;
    
    species->calcInverseDFTExp (out, fc);
    TS_ASSERT_VECTOR_APPROX (out, expected);
  }
  
  void testCalculateEIR () {
    // setUp() should have done all the hard work initialising stuff
    // current limitation: humans have no infectiousness
    //NOTE: could we load the model from a checkpoint at end of initialisation phase?
    
    // Run advancePeriod and getEIR for each human. Don't care about eir returned by getEIR;
    // just test vtm->eirPerDayOfYear and species parameters.
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (doc["calculateEIR"]["output"]);
  }
  
  void testCalculateEirDeterrency () {
    const YAML::Node& node = doc["calculateEirDeterrency"];
    species->ITNDeterrency = yaml2WeibullDecayedValue(node["Deterrency"]);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupITN();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirPreprandialKilling () {
    const YAML::Node& node = doc["calculateEirPreprandialKilling"];
    species->ITNPreprandialKillingEffect = yaml2WeibullDecayedValue(node["PreprandialKilling"]);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupITN();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirPostprandialKilling () {
    const YAML::Node& node = doc["calculateEirPostprandialKilling"];
    species->ITNPostprandialKillingEffect = yaml2WeibullDecayedValue(node["PostprandialKilling"]);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupITN();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirRestKilling () {
    const YAML::Node& node = doc["calculateEirRestKilling"];
    species->IRSKillingEffect = yaml2WeibullDecayedValue(node["RestKilling"]);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      it->setupIRS();
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void testCalculateEirLarviciding () {
    const YAML::Node& node = doc["calculateEirLarviciding"];
    node["Larviciding"]["effectiveness"] >> species->larvicidingIneffectiveness;
    species->larvicidingIneffectiveness = 1 - species->larvicidingIneffectiveness;
    species->larvicidingEndStep = 1000;
    vtm->timeStepNumEntoInnocs = 0;
    vtm->advancePeriod (*population, simulationTime);
    for (list<Human>::iterator it = population->begin(); it != population->end(); ++it)
      vtm->getEIR (simulationTime, it->perHostTransmission, it->getAgeInYears());
    TS_ASSERT_SPECIES_APPROX (node["output"]);
  }
  
  void doAssertSpecies (const char *f, unsigned l, const YAML::Node& node) {
    ETS_ASSERT_EQUALS ((int)vtm->timeStepNumEntoInnocs, simulation->_population->_populationSize);
    /* Print values so they can easily be copied into expected output:
    cout << endl << setprecision(10);
    cout << "average EIR: " << vtm->timeStepTotalEir / vtm->timeStepTotalEirEntries;
    cout << "\nP_A:\t" << species->P_A;
    cout << "\nP_df:\t" << species->P_df;
    cout << "\nP_dif:\t" << species->P_dif;
    cout << "\nN_v:\t" << species->N_v;
    cout << "\nO_v:\t" << species->O_v;
    cout << "\nS_v:\t" << species->S_v;
    */
    double avEIR;
    node["averageEIR"] >> avEIR;
    double result = 0.0;
    for (vector<double>::const_iterator it = vtm->timeStepEntoInnocs.begin(); it != vtm->timeStepEntoInnocs.end(); ++it)
      result += *it;
    TS_ASSERT_APPROX (result / vtm->timeStepNumEntoInnocs, avEIR);
    _TS_ASSERT_VECTOR_APPROX (f,l, species->P_A, yaml2Std<double> (node["P_A"]));
    _TS_ASSERT_VECTOR_APPROX (f,l, species->P_df, yaml2Std<double> (node["P_df"]));
    _TS_ASSERT_VECTOR_APPROX (f,l, species->P_dif, yaml2Std<double> (node["P_dif"]));
    _TS_ASSERT_VECTOR_APPROX (f,l, species->N_v, yaml2Std<double> (node["N_v"]));
    _TS_ASSERT_VECTOR_APPROX (f,l, species->O_v, yaml2Std<double> (node["O_v"]));
    _TS_ASSERT_VECTOR_APPROX (f,l, species->S_v, yaml2Std<double> (node["S_v"]));
  }
  
private:
  YAML::Node doc;
  int simulationTime;
  Simulation *simulation;
  VectorTransmission *vtm;
  VectorTransmissionSpecies *species;
  list<Human> *population;
};

#endif
