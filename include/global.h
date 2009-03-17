#ifndef Hmod_global
#define Hmod_global
#include "Constant.h"
#include <fcntl.h>
#include <math.h>
#include <vector>
using namespace std;

class Global
{
public:
  /// Sets parameters in Global.
  static void initGlobal ();
  
  /// Sets modelVersion, checking for incompatible versions.
  static void setModelVersion ();
  
  static int modIntervalsPerYear (int i);

  /// Variables that must be checkpointed.
  //@{

/** Model version defines which implementations of hard-coded options should be
 * used. The integer value of modelVersion passed from the .xml is converted to
 * binary with each bit corresponding to a different dichotomous option.  The
 * original default model is modelVersion=0 */
  static ModelVersion modelVersion;
  //@}
  
  /// Data read from xml which doesn't need to be checkpointed.
  //@{

  /// temporal resolution of simulation, in days
  static int interval;
   //Simulation time steps per year
  static int intervalsPerYear;
   //Maximum age of individuals in a scenario in time intervals
  static int maxAgeIntervals;
  static int simulationMode;
   //pre-erythrocytic latent period, in time steps
  static int latentp;
  
/*
  Size of the human population
  Moved from population.f so that transmission model.f can see it.
*/
  static vector<int> infantDeaths;
  static vector<int> infantIntervalsAtRisk;
  //@}
};

inline int isOptionIncluded (int allOptions, int option) {

  /*
    This is called soooo many timesthat performance optimization at the cost of readability is justified.
    isOptionIncluded=iand(shiftl(1,Option),AllOptions).gt.0
    return AllOptions AND ShiftLeft(1,Option Bits)
  */
  return allOptions & (1 << option);

};

#endif
