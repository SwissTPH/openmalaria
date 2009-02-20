/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#ifndef Hmod_human
#define Hmod_human
#include "drug.h"
#include "global.h"
#include "event.h"
#include <list>
#include "infection.h"
#include "withinHostModel.h"

extern double smuY;
//Standard dev innate immunity for densities
extern double sigma_i;
//Pyrogenic threshold at birth (Y*0)
extern double initPyroThres;
// Ystar2: critical value in determining increase in pyrogenic threshold
extern double Ystar2_13;
//alpha: factor determining increase in pyrogenic threshold
extern double alpha14;
//Ystar1: critical value of parasite density in determing increase in pyrog t
extern double Ystar1_26;
//sevMal: critical density for severe malaria episode (Y*B1)
extern double sevMal_21;
//Critical age for co-morbidity (for both severe and indirect)
extern double critAgeComorb_30;
//comorbidity prevalence at birth as a risk factor for severe
extern double comorbintercept_24;
//comorbidity prevalence at birth as a risk factor for indirect
extern double indirRiskCoFactor_18;
// contribution of parasite densities to acquired immunity in the presence of fever
extern double immPenalty_22;
/*
  Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
  This variable decays the effectors cumulativeH and cumulativeY in a way that their
  effects on densities (1-Dh and 1-Dy) decay exponentially.
*/
extern double asexImmRemain;
/*
  Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
  This variable decays the effectors cumulativeH and cumulativeY exponentially.
*/
extern double immEffectorRemain;
/*
  Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0
*/
extern double rateMultiplier_31;
extern double densityExponent_32;
extern double BaselineAvailabilityShapeParam;
/*
  The detection limit (in parasites/ul) is currently the same for PCR and for microscopy
  TODO: in fact the detection limit in Garki should be the same as the PCR detection limit
  The density bias allows the detection limit for microscopy to be higher for other sites
*/
extern double detectionlimit;

// Forward declaration
class CaseManagement;
class TransmissionModel;

//! Model of a human individual 
class Human {


 public:


  //! Constructor which does not need random numbers. Only for testing.
  Human();

  //!  Initialise all variables of a human datatype including infectionlist
  //!  and druglist
  /*  
      \param ID unique identifier
      \param dob date of birth in time steps
  */
  Human(int ID, int dateOfBirth, CaseManagement* caseManagement, int simulationTime);

  //!  Initialise all variables of a human datatype including infectionlist
  //! and druglist
  /*
    \param funit IO unit
  */
  Human(istream& funit, CaseManagement* caseManagement, int simulationTime);

  ~Human();

  friend ostream& operator<<(ostream& out, const Human& human);
  friend istream& operator>>(istream& in, Human& human);

  void update(int simulationTime, TransmissionModel* transmissionModel);

   //! Reads out from a file the state data of the current human
  /* 
     \param funit io unit number 
  */
  void readFromFile(fstream& funit);
  
  //! Clears all infections in an individual
  void clearAllInfections();

  /*! Apply interventions to this human if eligible. Calculate the remaining
      efficacy of the latest vaccination if vaccinated before */
  void updateInterventionStatus();

  /*! Until now, this only includes decay of immunity against
      asexual blood stages */
  void updateImmuneStatus();

  void updateTime(int simulationTime) { _simulationTime = simulationTime; };
  
  void updateInfection(double adultEIR, double availability, double expectedNumberOfInfections);

  void determineClinicalStatus();

  void treatInfections();

  void setWeight(double weight){ _weight=weight; };
  double getWeight(){ return _weight; };

  //! docu 
  void medicate(string drugName, double qty, int time);

  /*! Determine if this human experiences a clinical event, and if so, which
    type of event */
  void determineClinicalEvent();
  
  /*! 
    Update the number of doses and the date of the most recent vaccination in
    this human */
  void vaccinate();

  //! Summarize the state of a human individual.
  void summarize();

  //! Determines the age group of a human
  int ageGroup() const;

  //! Get the age in years, based on an input reference time.
  double getAgeInYears(int time) const;

  double getCumulativeYlag() const {return _cumulativeYlag;}
  void setCumulativeYlag(double lag) {_cumulativeYlag = lag;}

  //! Get the PEV Efficacy
  double getPEVEfficacy() {return _PEVEfficacy;} 

  //! Get the Availability to mosquitoes
  double getBaselineAvailabilityToMosquitoes() {return _BaselineAvailabilityToMosquitoes;};

  double getProbTransmissionToMosquito() {return _ptransmit;};

  //! Get the cumulative EIRa
  double getCumulativeEIRa() {return _cumulativeEIRa;};

  //! Return doomed value
  int getDoomed() {return _doomed;};

  //! Returns the date of birth
  int getDateOfBirth() {return _dateOfBirth;};

  double getProbabilityOfInfection() { return _pinfected;};

  //! Get the age in years, based on current simulationTime in human.
  double getAgeInYears() const;

  double getTimeStepMaxDensity() const {return _timeStepMaxDensity;}
  void setTimeStepMaxDensity(double timeStepMaxDensity){_timeStepMaxDensity = timeStepMaxDensity;}
  int getSimulationTime() const {return _simulationTime;}
  double getTotalDensity() const {return _totalDensity;}
  void setTotalDensity(double totalDensity) {_totalDensity=totalDensity;}
  double getCumulativeh() const {return _cumulativeh;}
  void setCumulativeh(double cumulativeh) {_cumulativeh = cumulativeh;}
  double getCumulativeY() const {return _cumulativeY;}
  void setCumulativeY(double cumulativeY) {_cumulativeY = cumulativeY;}
  void setPTransmit(double pTransmit) {_ptransmit = pTransmit;}
  void setPyrogenThres(double pyrogenThres) {_pyrogenThres = pyrogenThres;}
  double getInnateImmunity() const {return _innateImmunity;}
  int getPatentInfections() const {return _patentinfections;}
  void setPatentInfections(int patentinfections) {_patentinfections=patentinfections;}
  double getBSVEfficacy() {return _BSVEfficacy;}
  int getMOI() const {return _MOI;}
  void setMOI(int MOI) {_MOI=MOI;}
  std::list<Infection*>* getInfections() { return &infections;}
  int getLastIPTIorPlacebo() {return _lastIptiOrPlacebo;};
  void setLastIPTIorPlacebo(int last) {_lastIptiOrPlacebo = last;};
  int getLastSPDose() const {return _lastSPDose;};
  void setSPDose(int dose) {_lastSPDose=dose;};


  /*!  Determines the probability that the individual transmits to a feeding
    mosquito */
  double infectiousness();

  //! Returns Cumulative Infections
  int getCumulativeInfections() {return _cumulativeInfections;};

  //! Determine the current pyrogenic threshold.
  double Ystar();

  void setProbabilityOfInfection(double probability) { _pinfected=probability;};

  //! Return if human has an insecticide treated net
  bool hasInsecticideTreatedNet() {return _ITN;}

  //! Set doomed
  void setDoomed(int doomed) { _doomed = doomed;}

  //! Update the cumulative EIRa
  void updateCumulativeEIRa(double value) { _cumulativeEIRa+=value; };

  //! Set a new caseManagement. Used if the changeHS intervention is called.
  void setCaseManagement(CaseManagement* caseManagement);

  list<Drug*>* getDrugs();

  WithinHostModel* _withinHostModel;
 
 private:
  

  // Time from start of the simulation
  int _simulationTime;
  
  CaseManagement* _caseManagement;
  std::list<Infection*> infections;

  int _currentWHSlot;
  double _weight;
  //!Total asexual blood stage density
  double _ylag[4];
  //!last SP Dose given
  int _lastSPDose;
  //!last IPTi or placebo dose given 
  int _lastIptiOrPlacebo;
  //!Cumulative number of infections since birth
  int _cumulativeInfections;
  //!Date of birth, time step since start of warmup
  int _dateOfBirth;
  //!Indicates that individual will die from indirect mortality
  int _doomed;
  //!unique identifier
  int _ID;
  //!number of vaccine doses this individual has received
  int _lastVaccineDose;
  //!indicates the latest treatment regimen(1st, 2nd or 3rd line)
  int _latestRegimen;
  //!time of the last treatment
  int _tLastTreatment;
  //!multiplicity of infection
  int _MOI;
  //!Total asexual blood stage density
  double _totalDensity;
  //!Number of infections with densities above the limit of detection
  int _patentinfections;
  //!Remaining efficacy of Blood-stage vaccines
  double _BSVEfficacy;
  //!cumulativeY from previous timestep
  double _cumulativeYlag;
  //!Number of infective bites since birth
  double _cumulativeEIRa;
  //!Number of infections received since birth
  double _cumulativeh;
  //!Cumulative parasite density since birth
  double _cumulativeY;
  //!innate ability to control parasite densities
  double _innateImmunity;
  //!Remaining efficacy of Pre-erythrocytic vaccines
  double _PEVEfficacy;
  //!pinfected: probability of infection (cumulative or reset to zero in massTreatment). Appears to be used only for calculating expected inoculations for the analysis
  //!of pre-erythrocytic immunity.
  double _pinfected;
  //!probability that a mosquito will become infected if feeds on individual
  double _ptransmit;
  //!critical density for fever (clinical episodes)
  double _pyrogenThres;
  //!Remaining efficacy of Transmission-blocking vaccines
  double _TBVEfficacy;
  //!has an ITN (Insecticide treated net)
  bool _ITN;
  //!Maximum parasite density during the previous 5-day interval
  double _timeStepMaxDensity;
  //!Baseline availability to mosquitoes
  double _comorbidityFactor; 
  //! comorbidity factor for heterogeneity 
  double _treatmentSeekingFactor;
  //! treatment seeking for heterogeneity
  double _BaselineAvailabilityToMosquitoes;
  Event _latestEvent;
  list<Drug*> _drugs;
  DrugProxy* _proxy;

  //! Create a new infection requires that the human is allocated and current
  void newInfection();

  /*!  Clears all infections which have expired (their startdate+duration is less
      than the current time). */
  void clearOldInfections();

  //! Treats all infections in an individual
  void treatAllInfections();

  /*!  SP drug action applies to each infection depending on genotype and when
    the individual had their last dose of SP */
  void SPAction();

  //! Determines eligibility and gives IPTi SP or placebo doses 
  void setLastSPDose();

  void updateInfectionStatus();


   /*! should return true in case of effective or partially effective
     treatment, false otherwise */
  bool uncomplicatedEvent(bool isMalaria);

  /*!  Determines whether there is an acute episode, or concomitant fever and
    then whether the episode is severe, uncomplicated or there is an indirect
    death returns:true if there was an effective treatment, false otherwise.
  */
  bool defineEvent();
  
  //! returns true in case of effective treatment, false otherwise
  bool severeMalaria();
  
  //! Introduce new infections via a stochastic process
  void introduceInfections(double adultEIR, double availability, double expectedNumberOfInfections);

  //! docu 
  void computeFinalConc(Drug drg, int hl);

  //! In practice calculates drug concs (PK) for each drug
  void startWHCycle();

  double convertDose(Drug drg, Dose ds);


  //! docu
  double computeExponentialDecay(double c, int hl, int t);

  //! docu
  double calculateSelectionCoefficient(Infection inf);

  void clearInfection(Infection *iCurrent);

  /*!
    Demonstrates the use of the CM-related calls to medicate. May or may not be
    used as a subroutine in the final implementation.
  */
  void doCM(int entrypoint);

  //! Reads drugs from checkpoint
  void readDrugs(fstream& funit, int multiplicity);

  //!  Write all drugs of this human to the standard checkpoint file.
  /*
    \param funit io unit number
  */
  void writeDrugs(fstream& funit);


};

//Start DrugAction declarations
extern Human * currentHuman;
extern int currentWHSlot;
/*
  Current WH time slot
  End DrugAction declarations
*/
void initHumanParameters ();
void initDrugAction();

#endif


