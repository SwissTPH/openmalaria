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
#include "population.h"
#include "util/gsl.h"
#include "inputData.h"
#include "Human.h"
#include "simulation.h"
#include "intervention.h"
#include "Transmission/TransmissionModel.h"
#include "Clinical/ImmediateOutcomes.h"		// for changeHS intervention
#include "Transmission/NonVector.h"	// changeEIR intervention deals directly with this model
#include "summary.h"
#include "Pathogenesis/PathogenesisModel.h"
#include "NeonatalMortality.h"
#include <math.h>

using namespace std;


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

double *Population::cumpc;

int Population::IDCounter;

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
  Human::initHumanParameters();
  NeonatalMortality::init();
#ifdef OMP_CSV_REPORTING
  csvReporting.open ("population.csv", ios::app);
#endif
}

void Population::clear(){
  Human::clear();
#ifdef OMP_CSV_REPORTING
  csvReporting.close();
#endif
}


// -----  non-static methods: creation/destruction, checkpointing  -----

Population::Population()
    : _populationSize(get_populationsize())
{
  _transmissionModel = TransmissionModel::createTransmissionModel();

  _workUnitIdentifier=get_wu_id();
  _maxTimestepsPerLife=maxLifetimeDays/Global::interval;
  cumpc = (double*)malloc(((_maxTimestepsPerLife))*sizeof(double));

  IDCounter=0;
}

Population::~Population() {
  for(HumanIter iter=_population.begin(); iter != _population.end(); ++iter){
    iter->destroy();
  }
  delete _transmissionModel;  
}

void Population::read (istream& in) {
  //Start reading a checkpoint
  _transmissionModel->read (in);
  in >> _populationSize;
  if (_populationSize != get_populationsize())
    throw checkpoint_error("population size incorrect");
  in >> IDCounter;
  in >> mu0;
  in >> mu1;
  in >> alpha0;
  in >> alpha1;
  in >> rho;
  in >> _workUnitIdentifier;

  if (_workUnitIdentifier !=  get_wu_id()) {
    cerr << "cp_ct " << get_wu_id() << ", " << _workUnitIdentifier << endl;
    exit(-9);
  }

  //Start reading the human data
  int popSize;
  in >> popSize;
  if (popSize > _populationSize)
    throw checkpoint_error ("population size exceeds that given in scenario.xml");
  int indCounter = 0;	// Number of individuals read from checkpoint
  while(!(in.eof()||popSize==indCounter)){
    _population.push_back(Human(in, *_transmissionModel));
    indCounter++;
  }
  if (popSize != indCounter)
    throw checkpoint_error("can't read whole population (out of data)");
}
void Population::write (ostream& out) {
  _transmissionModel->write (out);
  out << _populationSize << endl;
  out << IDCounter << endl;
  out << mu0 << endl;
  out << mu1 << endl;
  out << alpha0 << endl;
  out << alpha1 << endl ;
  out << rho << endl; 
  out << _workUnitIdentifier << endl;

  //Write human data
  out << _population.size() << endl;	// this may not be equal to _populationSize due to startup optimisation
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    iter->write (out);
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
  const scnXml::AgeGroupPerC::GroupSequence& group = getDemography().getAgeGroup().getGroup();
  if (group.size() < ngroups-1) {
    ostringstream msg;
    msg << "expected " << ngroups-1 << " elements of \"group\" in demography->ageGroup (in scenario.xml)";
    throw xml_scenario_error(msg.str());
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
  rho = get_growthrate() * (0.01 * Global::yearsPerInterval);
  if (rho != 0.0)
    // Issue: in this case the total population size differs from _populationSize,
    // however, some code currently uses this as the total population size.
    throw xml_scenario_error ("Population growth rate provided.");
  
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
  // 1. Calculate cumpc
  cumpc[0] = 0.0;
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    double ageYears = (_maxTimestepsPerLife-j-1) * Global::yearsPerInterval;
    double M1s=(mu0 * (1.0-exp(-alpha0*ageYears)) / alpha0);
    double M2s=(mu1 * (exp(alpha1*ageYears)-1.0) / alpha1);
    double Ms=M1s+M2s;
    double predperc = exp(-rho*ageYears-Ms);
    if (j < _maxTimestepsPerLife-Global::maxAgeIntervals){
      predperc=0.0;
    }
    cumpc[j] = cumpc[j-1] + predperc;
  }
  int cumulativePop=0;
  double totalCumPC = cumpc[_maxTimestepsPerLife-1];
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    int iage=_maxTimestepsPerLife-j-1;
    //Scale using the total cumpc
    cumpc[j]=cumpc[j]/totalCumPC;
    // 2. Create humans
    if (!isCheckpoint){
      int targetPop = (int)floor(cumpc[j]*_populationSize+0.5);
      while (cumulativePop < targetPop) {
	if (InitPopOpt && iage > 0) {}	// only those with age 0 should be created here
	else newHuman(-iage);
	++cumulativePop;
      }
    }
  }
  
  // 3. Vector setup dependant on human population
  _transmissionModel->setupNv0 (_population, _populationSize);
}

void Population::preMainSimInit () {
  _transmissionModel->initMainSimulation();

  for (size_t i=0;i<Global::intervalsPerYear; i++) {
    ClinicalModel::infantIntervalsAtRisk[i]=0;
    ClinicalModel::infantDeaths[i]=0;
  }
}

// -----  non-static methods: simulation loop  -----

void Population::newHuman(int dob){
  ++IDCounter;
  _population.push_back(Human(*_transmissionModel, IDCounter, dob, Simulation::simulationTime));
}

void Population::update1(){
  NeonatalMortality::update (_population);
  // This should be called before humans contract new infections in the simulation step.
  _transmissionModel->vectorUpdate (_population, Simulation::simulationTime);
  
  //targetPop is the population size at time t allowing population growth
  int targetPop = (int)(_populationSize * exp(rho*Simulation::simulationTime));
  int cumPop = 0;
  
  // Update each human in turn
  //std::cout<<" time " <<t<<std::endl;
  HumanIter last = _population.end();
  --last;
  for (HumanIter iter = _population.begin(); iter != _population.end();){
    // Update human, and remove if too old:
    if (iter->update(Simulation::simulationTime,_transmissionModel)){
      iter->destroy();
      iter=_population.erase(iter);
      continue;
    }
    
    //BEGIN Population size & age structure
    ++cumPop;
    int age=(Simulation::simulationTime-iter->getDateOfBirth());
    
    // if (Actual number of people so far > target population size for this age) ...
    //FIXME: The +2 here is to replicate old results. I think it's wrong though. Also, it looks like this code assumes the maximum age of indivs is _maxTimestepsPerLife not Global::maxAgeIntervals.
    if (cumPop > targetCumPop (age+2, targetPop)) {
      --cumPop;
      iter->destroy();
      iter = _population.erase(iter);
      continue;
    }
    //END Population size & age structure
    ++iter;
  }	// end of per-human updates
  
  // increase population size to targetPop
  if (InitPopOpt && Simulation::simulationTime < Global::maxAgeIntervals) {
    // We only want people at oldest,
    //  Global::maxAgeIntervals - (Global::maxAgeIntervals-Simulation::simulationTime)
    // Hence pop size available is (total - popSize for anyone older):
    targetPop = targetPop - targetCumPop(Simulation::simulationTime+1, targetPop);
  }
  while (cumPop < targetPop) {
    newHuman(Simulation::simulationTime);
    //++nCounter;
    ++cumPop;
  }
  
  _transmissionModel->updateKappa (_population, Simulation::simulationTime);
  
#ifdef OMP_CSV_REPORTING
  if (Simulation::simulationTime % (Global::intervalsPerYear*5)==0) {
    csvReporting << Simulation::simulationTime << ',';
    list<Human>::reverse_iterator it = _population.rbegin();
    for (double ageLim = 0; ageLim <= maxLifetimeDays/365.0; ageLim += 1) {
      int counter=0;
      while (it != _population.rend() && it->getAgeInYears() < ageLim) {
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
  return (int)floor(cumpc[_maxTimestepsPerLife+1-ageTSteps] * targetPop + 0.5);
}


// -----  non-static methods: summarising and interventions  -----

void Population::newSurvey () {
  for(HumanIter iter=_population.begin(); iter != _population.end(); iter++){
    iter->summarize();
  }
  _transmissionModel->summarize (*Simulation::gMainSummary);
  Simulation::gMainSummary->incrementSurveyPeriod();
}

void Population::implementIntervention (int time) {
  const scnXml::Intervention* interv = getInterventionByTime (time);
  if (interv == NULL)
    return;
  
  if (interv->getChangeHS().present()) {
    if (Global::modelVersion & CLINICAL_EVENT_SCHEDULER)
      throw xml_scenario_error ("Only ClinicalImmediateOutcomes is compatible with change of health-system intervention.");
    changeHealthSystem (&interv->getChangeHS().get());
    ClinicalImmediateOutcomes::initParameters();	// should re-read all parameters
    
    //FIXME: surely we shouldn't do this at all? (DH)
    //TODO: Do we also need to re-init the kappa array?
    _transmissionModel->copyToInitialKappa();
  }
  
  if (interv->getChangeEIR().present()) {
    NonVectorTransmission* nvt = dynamic_cast<NonVectorTransmission*> (_transmissionModel);
    if (nvt != NULL) {
      nvt->setTransientEIR(interv->getChangeEIR().get());
    } else {
      throw xml_scenario_error("Warning: changeEIR intervention can only be used with NonVectorTransmission model!");
    }
  }
  
  if (interv->getVaccinate().present()) {
    massIntervention (interv->getVaccinate().get(), &Human::massVaccinate);
  }
  if (interv->getMDA().present()) {
    /* TODO: here we assume a 100% clearance rate for the MDA drug we use.
    This is not consistent with the way we treat according to the Health
    system description. The default clearance rate for MDA should be 100%
    since this simulates what was meant to happen in Garki.  We can change
    this by introducing an optional clearance rate that can be < 100%. */
    massIntervention(interv->getMDA().get(), &Human::clearInfections);
  }
  if (interv->getIpti().present()) {
    massIntervention (interv->getIpti().get(), &Human::IPTiTreatment);
  }
  
  if (interv->getITN().present()) {
    massIntervention (interv->getITN().get(), &Human::setupITN);
  }
  if (interv->getIRS().present()) {
    massIntervention (interv->getIRS().get(), &Human::setupIRS);
  }
  if (interv->getVectorAvailability().present()) {
    massIntervention (interv->getVectorAvailability().get(), &Human::setupVA);
  }
  
  if (interv->getLarviciding().present()) {
    _transmissionModel->intervLarviciding (interv->getLarviciding().get());
  }
}

void Population::massIntervention (const scnXml::Mass& mass, void (Human::*intervention)())
{
  double minAge = mass.getMinAge().present() ?
      mass.getMinAge().get() : 0.0;
  double maxAge = mass.getMaxAge().present() ?
      mass.getMaxAge().get() : 100.0;
  double coverage = mass.getCoverage();
  
  for(HumanIter iter=_population.begin(); iter != _population.end(); ++iter) {
    double ageYears = iter->getAgeInYears();
    if ((ageYears > minAge) && (ageYears < maxAge) && gsl::rngUniform() < coverage)
      // This is UGLY syntax. It just means call intervention() on the human pointed by iter.
      ((*iter).*intervention)();
  }
}
