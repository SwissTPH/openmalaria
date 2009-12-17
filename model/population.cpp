/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNu General Public License as published by
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
#include "Simulation.h"
#include "population.h"
#include "util/gsl.h"
#include "inputData.h"
#include "Surveys.h"

#include "Transmission/TransmissionModel.h"
#include "Transmission/NonVector.h"	// changeEIR intervention deals directly with this model

#include "Host/Human.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ImmediateOutcomes.h"		// for changeHS intervention
#include "Pathogenesis/PathogenesisModel.h"
#include "PkPd/PkPdModel.h"

#include "util/errors.hpp"
#include "util/ModelOptions.hpp"

#include <math.h>

namespace OM {
    
// -----  static data / methods  -----

double Population::ageGroupBounds[ngroups+1];
double Population::ageGroupPercent[ngroups];

double Population::M1[ngroups];
double Population::M2[ngroups];
double Population::M[ngroups];
double Population::pred[ngroups];

double Population::mu0;
double Population::mu1;
double Population::alpha0;
double Population::alpha1;
double Population::rho;
vector<double> Population::cumAgeProp;

int Population::IDCounter;
int Population::_maxTimestepsPerLife;
int Population::_workUnitIdentifier;

// Theoretical optimisation: don't include individuals who will die before the
// start of the main simulation in the initialisation.
// Two consequences:
// setRiskFromMaternalInfection uses some of these individuals. Could enforce
// this for initialisation, or just ignore (probably won't have much effect).
// And it's practially impossible to keep the same random number stream, so
// comparing results for testing isn't easy.
//
// No longer possible due to new Vector init.
const bool InitPopOpt = false;

#ifdef OMP_CSV_REPORTING
ofstream csvReporting;
#endif


void Population::init(){
  Host::Human::initHumanParameters();
  Host::NeonatalMortality::init();
  PkPd::PkPdModel::init();
#ifdef OMP_CSV_REPORTING
  csvReporting.open ("population.csv", ios::app);
#endif
  
  _workUnitIdentifier=InputData.get_wu_id();
  _maxTimestepsPerLife=maxLifetimeDays/Global::interval;
  cumAgeProp.resize (_maxTimestepsPerLife);
  
  IDCounter=0;
}

void Population::clear(){
    PkPd::PkPdModel::cleanup ();
  Host::Human::clear();
#ifdef OMP_CSV_REPORTING
  csvReporting.close();
#endif
}

void Population::staticRead (istream& in) {
  Host::NeonatalMortality::read (in);
  Clinical::ClinicalModel::staticRead(in);
  PkPd::PkPdModel::readStatic (in);
  
  in >> IDCounter;
  in >> mu0;
  in >> mu1;
  in >> alpha0;
  in >> alpha1;
  in >> rho;
}
void Population::staticWrite (ostream& out) {
  Host::NeonatalMortality::write (out);
  Clinical::ClinicalModel::staticWrite(out);
  PkPd::PkPdModel::writeStatic (out);
  
  out << IDCounter << endl;
  out << mu0 << endl;
  out << mu1 << endl;
  out << alpha0 << endl;
  out << alpha1 << endl ;
  out << rho << endl; 
}


// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population()
: populationSize(InputData.get_populationsize())
{
  _transmissionModel = Transmission::TransmissionModel::createTransmissionModel();
}

Population::~Population() {
  for(HumanIter iter=population.begin(); iter != population.end(); ++iter){
    iter->destroy();
  }
  delete _transmissionModel;  
}

void Population::checkpoint (istream& stream) {
    // Validate scenario.xml and checkpoint files correspond:
    _workUnitIdentifier & stream;
    if (_workUnitIdentifier !=  InputData.get_wu_id()) {
	cerr << "cp_ct " << InputData.get_wu_id() << ", " << _workUnitIdentifier << endl;
	exit(-9);
    }
    
    int popSize;
    popSize & stream;
    if (popSize != populationSize)
	throw util::checkpoint_error ("population size exceeds that given in scenario.xml");
    while(popSize > 0 && !stream.eof()){
	// Note: calling this constructor of Host::Human is slightly wasteful, but avoids the need for another
	// ctor and leaves less opportunity for uninitialized memory.
	population.push_back(Host::Human(*_transmissionModel, 0,0,0));
	population.back() & stream;
	--popSize;
    }
    if (int(population.size()) != populationSize)
	throw util::checkpoint_error("can't read whole population (out of data)");
}
void Population::checkpoint (ostream& stream) {
  _workUnitIdentifier & stream;

  //Write human data
  population.size() & stream;	// this may not be equal to populationSize due to startup optimisation
  for(HumanIter iter=population.begin(); iter != population.end(); ++iter){
    (*iter) & stream;
  }
}


// -----  static & non-static methods: further set up  -----

void Population::estimateRemovalRates () {
  // mu1, alpha1: These are estimated here
  // mu0, alpha0: These are fixed alpha0=4.0, mu0 calculated from the other parameters
  // rho - population growth rate (input)
  
  //double x[2];
  double rss;
  double tol;
  //double difference;
  double sumperc;
  //double w1[2];
  //double w2[2];
  //double w3[2];
  //double w4[2];
  //double w5[3];
  //double w6[3*2];
  int maxcal;
  int npar;
  int iw;
  double p1 = 0.371626412;
  double p2 = 0.841209593;

  //Get lower and upper age bounds for age groups and cumulative precentage of population from field data
  sumperc=0.0;
  const scnXml::AgeGroupPerC::GroupSequence& group = InputData.getDemography().getAgeGroup().getGroup();
  if (group.size() < ngroups-1) {
    ostringstream msg;
    msg << "expected " << ngroups-1 << " elements of \"group\" in demography->ageGroup (in scenario.xml)";
    throw util::xml_scenario_error(msg.str());
  }
  //Add age group for first month of life
  ageGroupBounds[0]=0.0;
  ageGroupBounds[1]=1.0/12.0;
  ageGroupPercent[0]=0.0;
  for (int i=1;i<ngroups; i++) {
    ageGroupBounds[i+1] = group[i-1].getUpperbound();
    ageGroupPercent[i] = group[i-1].getPoppercent();
    sumperc += ageGroupPercent[i];
  }
  sumperc = 100.0 / sumperc;	// multiplier to get percentages
  for (int i=0;i<ngroups; i++){
    ageGroupPercent[i]  = ageGroupPercent[i] * sumperc;
  }
    /*
  RSS between observed and predicted log percentage of population in age groups
  is minimised for values of mu1 and alpha1  
  calls setDemoParameters to calculate the RSS
    */
  tol=0.00000000001;
  maxcal=500000;
  npar=2;
  iw=3;
  rss=gsl::minimizeCalc_rss(&p1, &p2);
}

// Static method used by estimateRemovalRates
double Population::setDemoParameters (double param1, double param2) {
    rho = InputData.get_growthrate() * (0.01 * Global::yearsPerInterval);
  if (rho != 0.0)
    // Issue: in this case the total population size differs from populationSize,
    // however, some code currently uses this as the total population size.
    throw util::xml_scenario_error ("Population growth rate provided.");
  
  const double IMR=0.1;
  double M_inf=-log(1-IMR);
  
  mu1	= exp(param1) / 100;
  alpha1= exp(param2) / 100;
  alpha0= 4.0;
  mu0	= (M_inf - mu1 * (exp(alpha1*0.5) - 1) * alpha0) /
      (alpha1 * (1 - exp(-alpha0*0.5)));
  
  double sumpred=0.0;
  for (int i=0; i<ngroups-1; i++) {
    double midpt = (ageGroupBounds[i+1] + ageGroupBounds[i])*0.5;
    M1[i] = mu0 * (1.0 - exp(-alpha0*midpt)) / alpha0;
    M2[i] = mu1 * (exp(alpha1*midpt) - 1.0) / alpha1;
    M[i]  = M1[i] + M2[i];
    pred[i] = (ageGroupBounds[i+1] - ageGroupBounds[i])
        * exp(-rho*midpt - M[i]);
    sumpred += pred[i];
  }
  for (int i=0; i<ngroups-1; i++) {
    pred[i] = pred[i] / sumpred * 100.0;
  }
  double L_inf = exp(-rho*0.5 - M[1]);
  double M_nn  =-log(1.0 - 0.4*(1 - exp(-M[1])));
  double L1    = 1.0/12.0 * exp(-rho/24.0 - M_nn);
  double perc_inf = ageGroupPercent[0] + ageGroupPercent[1];
  ageGroupPercent[0] = perc_inf * L1 / L_inf;
  ageGroupPercent[1] = perc_inf - ageGroupPercent[0];
  
  double valsetDemoParameters=0.0;
  for (int i=0; i<ngroups-1; i++) {
    double residual = log(pred[i]) - log(ageGroupPercent[i]);
    valsetDemoParameters += residual*residual;
  }
  return valsetDemoParameters;
}

void Population::setupPyramid(bool isCheckpoint) {
  // 1. Calculate cumAgeProp
  cumAgeProp[0] = 0.0;
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    double ageYears = (_maxTimestepsPerLife-j-1) * Global::yearsPerInterval;
    double M1s=(mu0 * (1.0-exp(-alpha0*ageYears)) / alpha0);
    double M2s=(mu1 * (exp(alpha1*ageYears)-1.0) / alpha1);
    double Ms=M1s+M2s;
    double predperc = exp(-rho*ageYears-Ms);
    if (j < _maxTimestepsPerLife-Global::maxAgeIntervals){
      predperc=0.0;
    }
    cumAgeProp[j] = cumAgeProp[j-1] + predperc;
  }
  int cumulativePop=0;
  double totalCumPC = cumAgeProp[_maxTimestepsPerLife-1];
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    //Scale using the total cumAgeProp
    cumAgeProp[j]=cumAgeProp[j]/totalCumPC;
  }
  
  if (!isCheckpoint) {
    // 2. Create humans
    for (int j=1;j<_maxTimestepsPerLife; j++) {
      int iage=_maxTimestepsPerLife-j-1;
      int targetPop = (int)floor(cumAgeProp[j]*populationSize+0.5);
      while (cumulativePop < targetPop) {
	if (InitPopOpt && iage > 0) {}	// only those with age 0 should be created here
	else newHuman(-iage);
	++cumulativePop;
      }
    }
    
    // 3. Vector setup dependant on human population
    _transmissionModel->setupNv0 (population, populationSize);
  }
}

void Population::preMainSimInit () {
  _transmissionModel->initMainSimulation();

  for (size_t i=0;i<Global::intervalsPerYear; i++) {
    Clinical::ClinicalModel::infantIntervalsAtRisk[i]=0;
    Clinical::ClinicalModel::infantDeaths[i]=0;
  }
}

// -----  non-static methods: simulation loop  -----

void Population::newHuman(int dob){
  ++IDCounter;
  population.push_back(Host::Human(*_transmissionModel, IDCounter, dob, Global::simulationTime));
}

void Population::update1(){
  // Calculate relative availability correction, so calls from vectorUpdate,
  // etc., will have a mean of 1.0.
  double meanRelativeAvailability = 0.0;
  for (std::list<Host::Human>::const_iterator h = population.begin(); h != population.end(); ++h)
    meanRelativeAvailability += Transmission::PerHostTransmission::relativeAvailabilityAge (h->getAgeInYears());
  Transmission::PerHostTransmission::ageCorrectionFactor = populationSize / meanRelativeAvailability;
  
  Host::NeonatalMortality::update (population);
  // This should be called before humans contract new infections in the simulation step.
  _transmissionModel->vectorUpdate (population, Global::simulationTime);
  
  //targetPop is the population size at time t allowing population growth
  int targetPop = (int)(populationSize * exp(rho*Global::simulationTime));
  int cumPop = 0;
  
  // Update each human in turn
  //std::cout<<" time " <<t<<std::endl;
  HumanIter last = population.end();
  --last;
  for (HumanIter iter = population.begin(); iter != population.end();){
    // Update human, and remove if too old:
    if (iter->update(Global::simulationTime,_transmissionModel)){
      iter->destroy();
      iter=population.erase(iter);
      continue;
    }
    
    //BEGIN Population size & age structure
    ++cumPop;
    int age=(Global::simulationTime-iter->getDateOfBirth());
    
    // if (Actual number of people so far > target population size for this age) ...
    //FIXME: The +2 here is to replicate old results. I think it's wrong though. Also, it looks like this code assumes the maximum age of indivs is _maxTimestepsPerLife not Global::maxAgeIntervals.
    if (cumPop > targetCumPop (age+2, targetPop)) {
      --cumPop;
      iter->destroy();
      iter = population.erase(iter);
      continue;
    }
    //END Population size & age structure
    ++iter;
  }	// end of per-human updates
  
  // increase population size to targetPop
  if (InitPopOpt && Global::simulationTime < Global::maxAgeIntervals) {
    // We only want people at oldest,
    //  Global::maxAgeIntervals - (Global::maxAgeIntervals-Global::simulationTime)
    // Hence pop size available is (total - popSize for anyone older):
    targetPop = targetPop - targetCumPop(Global::simulationTime+1, targetPop);
  }
  while (cumPop < targetPop) {
    newHuman(Global::simulationTime);
    //++nCounter;
    ++cumPop;
  }
  
  _transmissionModel->updateKappa (population, Global::simulationTime);
  
#ifdef OMP_CSV_REPORTING
  if (Global::simulationTime % (Global::intervalsPerYear*5)==0) {
    csvReporting << Global::simulationTime << ',';
    list<Host::Human>::reverse_iterator it = population.rbegin();
    for (double ageLim = 0; ageLim <= maxLifetimeDays/365.0; ageLim += 1) {
      int counter=0;
      while (it != population.rend() && it->getAgeInYears() < ageLim) {
	++counter;
	++it;
      }
      csvReporting << counter << ',';
    }
    csvReporting << endl;
  }
#endif
}

int Population::targetCumPop (int ageTSteps, int targetPop) {
  return (int)floor(cumAgeProp[_maxTimestepsPerLife+1-ageTSteps] * targetPop + 0.5);
}


// -----  non-static methods: summarising and interventions  -----

void Population::newSurvey () {
  Survey& current = *Surveys.current;
  for(HumanIter iter=population.begin(); iter != population.end(); iter++){
    iter->summarize(current);
  }
  _transmissionModel->summarize (current);
}

void Population::implementIntervention (int time) {
    const scnXml::Intervention* interv = InputData.getInterventionByTime (time);
  if (interv == NULL)
    return;
  
  // Given an intervention descriptor for this time point, check which
  // interventions are included. All should check when data loading that they
  // have data if used according to getActiveInterventions().
  
  if (interv->getChangeHS().present()) {
    if (util::ModelOptions::option (util::CLINICAL_EVENT_SCHEDULER))
      throw util::xml_scenario_error ("Only ClinicalImmediateOutcomes is compatible with change of health-system intervention.");
    InputData.changeHealthSystem (&interv->getChangeHS().get());
    Clinical::ClinicalImmediateOutcomes::initParameters();	// should re-read all parameters
    
    //FIXME: surely we shouldn't do this at all? (DH)
    //TODO: Do we also need to re-init the kappa array?
    _transmissionModel->copyToInitialKappa();
  }
  
  if (interv->getChangeEIR().present()) {
    _transmissionModel->changeEIRIntervention (interv->getChangeEIR().get());
  }
  
  if (interv->getVaccinate().present()) {
    massIntervention (interv->getVaccinate().get(), &Host::Human::massVaccinate);
  }
  if (interv->getMDA().present()) {
    /* TODO: here we assume a 100% clearance rate for the MDA drug we use.
    This is not consistent with the way we treat according to the Health
    system description. The default clearance rate for MDA should be 100%
    since this simulates what was meant to happen in Garki.  We can change
    this by introducing an optional clearance rate that can be < 100%. */
    massIntervention(interv->getMDA().get(), &Host::Human::clearInfections);
  }
  if (interv->getIpti().present()) {
    massIntervention (interv->getIpti().get(), &Host::Human::IPTiTreatment);
  }
  
  if (interv->getITN().present()) {
    massIntervention (interv->getITN().get(), &Host::Human::setupITN);
  }
  if (interv->getIRS().present()) {
    massIntervention (interv->getIRS().get(), &Host::Human::setupIRS);
  }
  if (interv->getVectorAvailability().present()) {
    massIntervention (interv->getVectorAvailability().get(), &Host::Human::setupVA);
  }
  
  if (interv->getLarviciding().present()) {
    _transmissionModel->intervLarviciding (interv->getLarviciding().get());
  }
}

void Population::massIntervention (const scnXml::Mass& mass, void (Host::Human::*intervention)())
{
  double minAge = mass.getMinAge().present() ?
      mass.getMinAge().get() : 0.0;
  double maxAge = mass.getMaxAge().present() ?
      mass.getMaxAge().get() : 100.0;
  double coverage = mass.getCoverage();
  
  for(HumanIter iter=population.begin(); iter != population.end(); ++iter) {
    double ageYears = iter->getAgeInYears();
    if ((ageYears > minAge) && (ageYears < maxAge) && gsl::rngUniform() < coverage)
      // This is UGLY syntax. It just means call intervention() on the human pointed by iter.
      ((*iter).*intervention)();
  }
}
}
