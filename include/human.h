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
#include "global.h"
#include "event.h"
#include "Infection.h"
#include "withinHostModel.h"
#include "EntoIntervention.h"
#include "morbidityModel.h"
#include "drug.h"

#include <list>

// Forward declaration
class CaseManagementModel;
class TransmissionModel;

//! Model of a human individual 
class Human {
public:
  
  /// Constructors
  //@{
  //! Constructor which does not need random numbers. Only for testing.
  Human();

  /** Initialise all variables of a human datatype including infectionlist and
   * druglist.
   * 
   * \param ID unique identifier
   * \param dateOfBirth date of birth in time steps */
  Human(int ID, int dateOfBirth, CaseManagementModel* caseManagement, int simulationTime);

  /**  Initialise all variables of a human datatype including infectionlist and
   * and druglist.
   * \param funit IO unit */
  Human(istream& funit, CaseManagementModel* caseManagement, int simulationTime);
  //@}
  
  /** Destructor
   * 
   * NOTE: this destructor does nothing to allow shallow copying to the
   * population list. Human::destroy() does the real freeing and must be
   * called explicitly. */
  ~Human() {}
  
  /// The real destructor
  void destroy();
  
  /// Checkpointing functions
  //@{
  friend ostream& operator<<(ostream& out, const Human& human);
  friend istream& operator>>(istream& in, Human& human);

  /** Reads out from a file the state data of the current human.
   * \param funit io unit number */
  void readFromFile(fstream& funit);
  //@}
  
  /// Functions not called from within Human but calling functions in Human
  //@{
  /// Updates an individual for a time-step.
  void update(int simulationTime, TransmissionModel* transmissionModel);
  
  //! Summarize the state of a human individual.
  void summarize();
  //@}
  
  /// updateInfection related functions
  //@{
  void updateInfection(double expectedInfectionRate, double expectedNumberOfInfections);
  //@}
  
  /// treatInfections related functions
  //@{
  void treatInfections();
  //@}
  
  /// updateImmuneStatus related functions
  //@{
  /*! Until now, this only includes decay of immunity against
      asexual blood stages */
  void updateImmuneStatus();
  //@}
  
  /// determineClinicalStatus related functions
  //@{
  void determineClinicalStatus();
  
  //! Clears all infections in an individual
  void clearAllInfections();
  
  //! docu 
  void medicate(string drugName, double qty, int time);
  //@}
  
  /// updateInterventionStatus related functions
  //@{
  /*! Apply interventions to this human if eligible. Calculate the remaining
      efficacy of the latest vaccination if vaccinated before */
  void updateInterventionStatus();
  
  /*! 
    Update the number of doses and the date of the most recent vaccination in
    this human */
  void vaccinate();
  //@}
  
  /// Functions used internally by more than one category above
  /// (excluding update and summarize)
  //@{
  //! Determines the age group of a human
  int ageGroup() const;
  
  //! Get the age in years, based on an input reference time.
  double getAgeInYears(int time) const;

  //! Get the age in years, based on current simulationTime in human.
  double getAgeInYears() const;
  //@}
  
  /// Other functions not called within Human
  //@{
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
  
  double getTotalDensity() const {return _totalDensity;}
  
  double getCumulativeY() const {return _cumulativeY;}
  
  void setLastIPTIorPlacebo(int last) {_lastIptiOrPlacebo = last;};
  
  void setSPDose(int dose) {_lastSPDose=dose;};

  //! Returns Cumulative Infections
  int getCumulativeInfections() {return _cumulativeInfections;};

  void setProbabilityOfInfection(double probability) { _pinfected=probability;};

  //! Set doomed
  void setDoomed(int doomed) { _doomed = doomed;}

  //! Set a new caseManagement. Used if the changeHS intervention is called.
  void setCaseManagement(CaseManagementModel* caseManagement);

  /** Availability of host to mosquitoes (α_i). */
  double entoAvailability () const;
  /** Probability of a mosquito succesfully biting a host (P_B_i). */
  double probMosqSurvivalBiting () const;
  /** Probability of a mosquito succesfully finding a resting
   * place and resting (P_C_i * P_D_i). */
  double probMosqSurvivalResting () const;
  //@}
  
  /// Functions only used by oldWithinHostModel.cpp
  //@{
  double getTimeStepMaxDensity() const {return _timeStepMaxDensity;}
  void setTimeStepMaxDensity(double timeStepMaxDensity){_timeStepMaxDensity = timeStepMaxDensity;}
  
  void setTotalDensity(double totalDensity) {_totalDensity=totalDensity;}
  
  double getCumulativeh() const {return _cumulativeh;}
  void setCumulativeh(double cumulativeh) {_cumulativeh = cumulativeh;}
  
  void setCumulativeY(double cumulativeY) {_cumulativeY = cumulativeY;}
  
  void setPTransmit(double pTransmit) {_ptransmit = pTransmit;}
  
  double getInnateImmunity() const {return _innateImmunity;}
  
  int getPatentInfections() const {return _patentinfections;}
  void setPatentInfections(int patentinfections) {_patentinfections=patentinfections;}
  
  double getBSVEfficacy() {return _BSVEfficacy;}
  
  int getMOI() const {return _MOI;}
  void setMOI(int MOI) {_MOI=MOI;}
  
  std::list<Infection*>* getInfections() { return &infections;}
  
  int getLastSPDose() const {return _lastSPDose;};
  
  /*!  Determines the probability that the individual transmits to a feeding
    mosquito */
  double infectiousness();

  //@}
  
  /// static public
  //@{
  static void initHumanParameters ();
  
/*
  The detection limit (in parasites/ul) is currently the same for PCR and for microscopy
  TODO: in fact the detection limit in Garki should be the same as the PCR detection limit
  The density bias allows the detection limit for microscopy to be higher for other sites
*/
  static double detectionlimit;
  
  /* Shape constant of (Gamma) distribution of availability
  real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
  static double BaselineAvailabilityShapeParam;
  
  /* nwtgrps is the number of age groups for which expected weights are read in
  for use in the age adjustment of the EIR.
  It is used both by Human and TransmissionModel. */
  static const int nwtgrps= 27; 
  //@}

private:
  ///Private variables
  //@{
  WithinHostModel* _withinHostModel;

  // Time from start of the simulation
  int _simulationTime;
  
  CaseManagementModel* _caseManagement;
  std::list<Infection*> infections;

  MorbidityModel* _morbidityModel;

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
  double _totalDensity;		// possibly move to WithinHostModel
  //!Number of infections with densities above the limit of detection
  int _patentinfections;	//TODO: move to WithinHostModel
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
  //!Remaining efficacy of Transmission-blocking vaccines
  double _TBVEfficacy;
  //!Maximum parasite density during the previous 5-day interval
  double _timeStepMaxDensity;	// WithinHostModel, used by Morbidity
  //!Baseline availability to mosquitoes
  double _comorbidityFactor; 
  //! comorbidity factor for heterogeneity 
  double _treatmentSeekingFactor;
  //! treatment seeking for heterogeneity
  double _BaselineAvailabilityToMosquitoes;
  Event _latestEvent;
  
  /// Encapsulates drug code for each human
  DrugProxy _proxy;
  
  /// Rate/probabilities before interventions. See functions.
  double _entoAvailability;
  double _probMosqSurvivalBiting;
  double _probMosqSurvivalResting;
  
  /// Intervention: an ITN (active if netEffectiveness > 0)
  EntoInterventionITN entoInterventionITN;
  /// Intervention: IRS (active if insecticide != 0)
  EntoInterventionIRS entoInterventionIRS;
  //@}
  
  /// updateInfection related functions
  //@{
  //! Create a new infection requires that the human is allocated and current
  void newInfection();

  /*!  Clears all infections which have expired (their startdate+duration is less
      than the current time). */
  void clearOldInfections();
  //@}
  
  /// treatInfections related functions
  //@{
  //! Treats all infections in an individual
  void treatAllInfections();
  //@}

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
  void introduceInfections(double expectedInfectionRate, double expectedNumberOfInfections);

  //! docu 
  void computeFinalConc(Drug drg, int hl);

  //! In practice calculates drug concs (PK) for each drug
  void startWHCycle();

  double convertDose(Drug drg, Dose ds);


  //! docu
  double computeExponentialDecay(double c, int hl, int t);

  //! docu
  double calculateSelectionCoefficient(Infection& inf);

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


  /* Static private */
  
//Standard dev innate immunity for densities
  static double sigma_i;
//sevMal: critical density for severe malaria episode (Y*B1)
  static double sevMal_21;
//Critical age for co-morbidity (for both severe and indirect)
  static double critAgeComorb_30;
//comorbidity prevalence at birth as a risk factor for severe
  static double comorbintercept_24;
//comorbidity prevalence at birth as a risk factor for indirect
  static double indirRiskCoFactor_18;
// contribution of parasite densities to acquired immunity in the presence of fever
  static double immPenalty_22;
/*
  Remaining immunity against asexual parasites(after time step, each of 2 components y and h)
  This variable decays the effectors cumulativeH and cumulativeY in a way that their
  effects on densities (1-Dh and 1-Dy) decay exponentially.
*/
  static double asexImmRemain;
/*
  Remaining immunity against asexual parasites(after each time step, each of 2 components y and h)
  This variable decays the effectors cumulativeH and cumulativeY exponentially.
*/
  static double immEffectorRemain;
  
  //! Relative weights by age group
  /** Relative weights, based on data in InputTables\wt_bites.csv 
  The data are for Kilombero, Tanzania, taken from the Keiser et al (diploma
  thesis). The original source was anthropometric studies by Inez Azevedo Reads
  in weights by age group. The weights are expressed as proportions of 0.5*those
  in the reference age group. */
  static const double wtprop[nwtgrps];
  
  static double baseEntoAvailability;
  static double baseProbMosqSurvivalBiting;
  static double baseProbMosqSurvivalResting;
};

#endif
