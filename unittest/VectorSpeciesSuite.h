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
#include <cxxtest/TestSuite.h>
#include "ExtraAsserts.h"
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "Transmission/VectorSpecies.h"


template<class T>
void operator>> (const YAML::Node& node, vector<T>& vec) {
  /* doesn't work! 
  vec.resize (node.size());
  for (size_t i = 0; i < node.size(); ++i) {
    cout << "iteartion: " << i<<endl;
    node[i] >> vec[i];
  }*/
  vec.resize (node.size());
  vec.resize (0);
  T t;
  for (YAML::Iterator it = node.begin(); it != node.end(); ++it) {
    *it >> t;
    vec.push_back (t);
  }
}


/** Unit test for VectorTransmissionSpecies code.
 *
 * Possible tests to write:
 * 
VectorSpecies unit testing:
  initialise:
    read parameters from XML
    calculate EIR array from fourier series, rotate and add to return param
    read emergence params from XML and scale
  destroy: (does nothing)
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
  calcInverseDFTExp:
    calculates an array from a fourier series
  rotateArray:
    offsets all element indicies by some x mod length
 */
class VectorSpeciesSuite : public CxxTest::TestSuite
{
public:
  VectorSpeciesSuite() {
    ifstream file("VectorSpeciesSuite.yaml");
    
    YAML::Parser parser(file);
    parser.GetNextDocument(doc);
    
    file.close();
  }
  
  /** Run first to check data is read OK.
   *
   * Note that using doc["key"] to access map elements isn't speed-optimal when
   * many key-value pairs exist, but for a small number it should be fine. */
  void testSuiteData() {
    string suite;
    doc["suite"] >> suite;
    TS_ASSERT_EQUALS (suite, "VectorSpeciesSuite");
  }
  void testCalcInverseDFTExp () {
    const YAML::Node& node = doc["calcInverseDFTExp"];
    vector<double> out(10), fc, expected;
    node["fourierCoefficients"] >> fc;
    node["output"] >> expected;
    
    species.calcInverseDFTExp (out, fc);
    TS_ASSERT_VECTOR_APPROX (out, expected);
  }
  
  /* Use to print out data.
  struct PrintStruct {
    void operator() (double x) {
      cout << x << ", ";
    }
  } printFunctor;
    cout << endl << setprecision(10);
    for_each (out.begin(), out.end(), printFunctor);
  */
  
private:
  YAML::Node doc;
  VectorTransmissionSpecies species;
};

#endif
