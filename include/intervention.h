#ifndef Hmod_intervention
#define Hmod_intervention
#include <fcntl.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
const int no_intervention= 0;
const int irs_intervention= 1;
const int mda_intervention= 2;
const int vaccine_intervention= 3;
const int change_eir_intervention= 4;
const int change_hs_intervention= 5;
const int ipti_intervention= 6;
const int preerythrocytic_reduces_h= 1;
const int erythrocytic_reduces_y= 2;
const int transmission_blocking_reduces_k= 3;

/*
  Vaccine specific parameters
  The vaccine type, as binary encoded number
*/
extern int vaccineType;
 

   /*
Target age for EPI-like vaccination, in time steps
TODO: Should be an integer?
*/
extern double *targetagetstep;
extern int targetagetstepX;
   //Coverage , as a proportion of the poulation in the target age range
extern double *vaccineCoverage;
extern int vaccineCoverageX;

/*
Vaccine type specific parameters
Initial mean efficacy, definition depends on vaccine type
*/
extern double *PEVInitialMeanEfficacy;
extern int PEVInitialMeanEfficacyX;
   extern double *BSVInitialMeanEfficacy;
extern int BSVInitialMeanEfficacyX;
   extern double *TBVInitialMeanEfficacy;
extern int TBVInitialMeanEfficacyX;
   //Distribution of efficacies among individuals, parameter to sample from beta dist.
extern double PEVefficacyB;
   extern double BSVefficacyB;
   extern double TBVefficacyB;
   //Decay rate
extern double PEVdecay;
   extern double BSVdecay;
   extern double TBVdecay;
   /*

ITN specific parameters
ITNs present or not
*/
extern short ITN;

//decay of ITNs
extern double ITNdecay;
extern double Pu0;
extern double Pu1;
extern double c;
extern double z;

/*
IPT specific parameters
IPT present or not
*/
extern short IPT;
   //Number of IPTi doses
extern int numberOfIPTiDoses;
   //Target age for IPTi doses, in time steps
extern int *iptiTargetagetstep;
extern int iptiTargetagetstepX;
   //Coverage , as a proportion of the poulation in the target age range
extern double *iptiCoverage;
extern int iptiCoverageX;
   //Values   
extern double iptiEffect;
   /*
These are actually not IPTi related values but related
to the genotypes: Since the genotype variability was first studied in a IPTi intervention,
they are defined here. (Logically, they belong to infection.f)
*/
extern int numberOfGenoTypes;
   extern double *genotypeFreq;
extern int genotypeFreqX;
   extern double *genotypeACR;
extern int genotypeACRX;
   extern int *genotypeProph;
extern int genotypeProphX;
   extern int *genotypeTolPeriod;
extern int genotypeTolPeriodX;
   extern double *genotypeAtten;
extern int genotypeAttenX;

  /*!
      Commen to all vaccine types. Number of vaccine doses that are given either
      through EPI or as EPI Boosters
   */
 extern int _numberOfEpiDoses;

  /*!
    Number of initial efficacy values in the XML. If more than this number or
    doses are given, we assume the vaccine efficacy goes to the last defined
    initial efficacy.
  */
extern int _numberOfInitEff;

// Initialization routines
void initInterventionParameters ();
  void initVaccineParameters ();
  void initITNParameters ();
  void initIPTIParameters ();
  
  // Destruction routines
  void clearVaccineParameters ();
  void clearITNParameters ();
  void clearIPTIParameters ();
void clearInterventionParameters ();


class Intervention{

 public:

  //end IPTi things
  Intervention();
  ~Intervention();

 private:
 

   /*!
      Commen to all vaccine types. Number of vaccine doses that are given either
      through EPI or as EPI Boosters
   */
  int _numberOfEpiDoses;

  /*!
    Number of initial efficacy values in the XML. If more than this number or
    doses are given, we assume the vaccine efficacy goes to the last defined
    initial efficacy.
  */
  int _numberOfInitEff;

};

#endif
