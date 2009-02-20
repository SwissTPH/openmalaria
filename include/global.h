#ifndef Hmod_global
#define Hmod_global
#include <fcntl.h>
#include <math.h>
#include <vector>
using namespace std;

namespace Diagnosis {

  enum Value { NON_MALARIA_FEVER,
               UNCOMPLICATED_MALARIA,
               SEVERE_MALARIA,
               INDIRECT_MALARIA_DEATH };
}

namespace Outcome {

/*
  Possibilities for outcomes are:
  for non-treated
*/
  enum Value {
    
    // non treated
    NO_CHANGE_IN_PARASITOLOGICAL_STATUS_NON_TREATED,
    //for outpatients
    NO_CHANGE_IN_PARASITOLOGICAL_STATUS_OUTPATIENTS,
    //for inpatients
    NO_CHANGE_IN_PARASITOLOGICAL_STATUS_INPATIENTS,
    //for non-treated
    PARASITES_ARE_CLEARED_PATIENT_RECOVERS_NON_TREATED,
    //for outpatients
    PARASITES_ARE_CLEARED_PATIENT_RECOVERS_OUTPATIENTS,
    //for inpatients
    PARASITES_ARE_CLEARED_PATIENT_RECOVERS_INPATIENTS,
    //for non-treated
    PARASITES_ARE_CLEARED_PATIENT_HAS_SEQUELAE_NON_TREATED,
    //for inpatients
    PARASITES_ARE_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS,
    //for non-treated
    PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_NON_TREATED,
    //for inpatients
    PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS,
    //for non-treated
    PATIENT_DIES_NON_TREATED,
    //for inpatients
    PATIENT_DIES_INPATIENTS,
    INDIRECT_DEATH,
    //for outpatients in models of pk/pD
    PARASITES_PKPD_DEPENDENT_RECOVERS_OUTPATIENTS
  };

}

const int  missing_value= -99999;
   //Days in a year
const int daysInYear= 365;
   /*
   TODO: Add simulation mode for entomological model.
   This probably needs to be done soon.
There are 3 simulation modes.
Equilibrium mode 
This is used for the warm-up period and if we want to separate direct effect of an intervention
from indirect effects via transmission intensity. The seasonal pattern and intensity of the EIR 
do not change over years 
*/
const int equilibriumMode= 2;
   /*
Transient EIR known
This is used to simulate an intervention that changes EIR, and where we have measurements of the EIR over 
time during the intervention period.
*/
const int transientEIRknown= 3;
   /*
EIRchanges
The EIR changes dynamically during the intervention phase as a function of the characteristics of the interventions
*/
const int dynamicEIR= 4;
   //The mean of the base line availability, is used by both human.f and entomology.f
const double BaselineAvailabilityMean= 1.0;
   /*
The relative risk of non malaria fever is used in human.f
TODO: This is currently set arbitrarily to 1.0, but should be defined in the .xml
*/
const double RelativeRiskNonMalariaFever= 1.0;
   /*

All variables that must be checkpointed. These are all values of which the variables 
 we need to resume state.
Model version defines which implementations of hard-coded options should be used
The integer value of modelVersion passed from the .xml is converted to binary
with each bit corresponding to a different dichotomous option.  The original default
model is modelVersion=0  
*/
extern int modelVersion;
   /*
TODO: can we document under each modelVersion what that model version does?
TODO: can we implement error trapping of incompatible model versions?
*/
const int penalisationEpisodes= 1;
   /*
Effective cumulative exposure to blood stage parasites is reduced during a clinical
episode, so that clinical episodes have a negative effect on blood stage immunity: 
default is no effect  
*/
const int negativeBinomialMassAction= 2;
   /*
Baseline availability of humans is sampled from a gamma distribution: 
default is to have fixed availability, and a decreasing proportion of inocula
surviving as EIR increases
*/
const int attenuationAsexualDensity= 3;
const int lognormalMassAction= 4;
   /*
Baseline availability of humans is sampled from a log normal distribution: 
default is to have fixed availability, and a decreasing proportion of inocula
surviving as EIR increases
*/
const int lognormalMassActionPlusPreImm= 5;
const int  maxDensCorrection= 6;
   /*
TODO: what is this exactly?
Possibly a better model version ensuring that the effect of variation in innate immunity
is reflected in case incidence:
*/
const int innateMaxDens= 7;
const int maxDensReset= 8;
   /*
Parasite density is determined by a dynamic intra host model:
default is to use empirical summary of malariatherapy densities
*/
const int withinHostParasite= 9;
   /*
Clinical episodes occur if parasitaemia exceeds the pyrogenic threshold: 
Default is to trigger clinical episodes stochastically, with the pyrogenic threshold
determining the density at which the probability is 50%
*/
const int predeterminedEpisodes= 10;
   /*
Simulates non-malaria fever attacks
Default is not to simulate non-malaria fever attacks
*/
const int nonMalariaFevers= 11;
   /*
Simulates PK and PD of drugs
Default is to treat drug action as binary
*/
const int includesPKPD= 12;
   /*
Uses version 2 of the case management model
Default is to use the Tediosi et al case management model
*/
const int caseManagementV2= 13;
 //Default is to use the Ross et al morbidity model
const int MuellerMorbidityModel= 14;
 //Allow simple heterogeneity in transmission, comorbidity and treatment-seeking
const int transHet = 15;  
const int comorbHet = 16;
const int treatHet = 17;
const int comorbTransHet = 18;
const int transTreatHet = 19;
const int comorbTreatHet = 20;
const int tripleHet = 21;
const int noVectorControl = 22; 

/*

Variables that need not be checkpointed. Note that when a simulation is resumed, the readInterface
subroutine is called, so input data does not have to be checkpointed.

temporal resolution of simulation, in days
*/
extern int interval;
   //Simulation time steps per year
extern int intervalsPerYear;
   //Maximum age of individuals in a scenario in time intervals
extern int maxAgeIntervals;
   extern int simulationMode;
   //pre-erythrocytic latent period, in time steps
extern double latentp;
  
/*
   Size of the human population
   Moved from population.f so that transmission model.f can see it.
*/
extern vector<int> infantDeaths;
extern int infantDeathsSize;
extern vector<int> infantIntervalsAtRisk;
extern int infantIntervalsAtRiskSize;


void initGlobal ();
void clearGlobalParameters ();
int modIntervalsPerYear (int i);

inline int isOptionIncluded (int allOptions, int option) {

  /*
    This is called soooo many timesthat performance optimization at the cost of readability is justified.
    isOptionIncluded=iand(shiftl(1,Option),AllOptions).gt.0
    return AllOptions AND ShiftLeft(1,Option Bits)
  */
  return allOptions & (1 << option);

};

double mymodf(double, double);


#define W_MINIMIZE_CALC_RSS w_minimize_calc_rss_


#ifdef _WIN32
int nearbyint(double x);
int round(double x);
#endif

#endif

