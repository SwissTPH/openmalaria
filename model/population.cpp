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
#include "GSLWrapper.h"
#include "timer.h"
#include "inputData.h"
#include "human.h"
#include "simulation.h"
#include "intervention.h"
#include "TransmissionModel.h"
#include "TransmissionModel/NonVector.h"	// changeEIR intervention deals directly with this model
#include "summary.h"
#include "Pathogenesis/PathogenesisModel.h"
#include "CaseManagementModel.h"
#include "BoincWrapper.h"
#include <math.h>

using namespace std;


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

Population::Population()
    : _populationSize(get_populationsize())
{
  _transmissionModel = TransmissionModel::createTransmissionModel();

  CaseManagementModel::init();
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

void Population::init(){
  Human::initHumanParameters();
}

void Population::clear(){
  Human::clear();
}

void Population::preMainSimInit () {
  _transmissionModel->initMainSimulation(_populationSize);

  initialiseInfantArrays();
}

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
  rss=w_minimize_calc_rss(&p1, &p2);
}

void Population::setupPyramid(bool isCheckpoint){
  cumpc[0] = 0.0;
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    double ageYears = (_maxTimestepsPerLife-j-1) / (1.0*Global::intervalsPerYear);
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
    if (!isCheckpoint){
      int N_new=(int)floor(cumpc[j]*_populationSize+0.5)-cumulativePop;
      for (int j1=1;j1<=N_new; j1++){
	newHuman(-iage);
	++cumulativePop;
      }
    }
  }
}

void Population::initialiseInfantArrays(){
  for (size_t i=0;i<Global::intervalsPerYear; i++) {
    Global::infantIntervalsAtRisk[i]=0;
    Global::infantDeaths[i]=0;
  }
}

void Population::initialiseHumanList(){
  IDCounter=0;
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
  start_cp_timer();
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    out << *iter;
  }

  //Finished writing lists
  stop_cp_timer();
}

void Population::read (istream& in) {
  //Start reading a checkpoint
  initialiseHumanList();
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
  int indCounter = 0;	// Number of individuals read from checkpoint
  while(!(in.eof()||_populationSize==indCounter)){
      //continue: Fortran cont is probably not C cont
    _population.push_back(Human(in, *_transmissionModel));
    indCounter++;
  }
  if (_populationSize != indCounter)
    throw checkpoint_error("can't read whole population (out of data)");
}

void Population::newHuman(int dob){
  ++IDCounter;
  _population.push_back(Human(*_transmissionModel, IDCounter, dob, Simulation::simulationTime));
}

void Population::update1(){
  //NOTE: this needs to be called somewhere; and should really be called before humans contract new infections in the simulation step
  //_transmissionModel->advancePeriod (_population, Simulation::simulationTime);
  
  int nCounter=0;	//NCounter is the number of indivs per demogr age group
  int pCounter=0;	//PCounter is the number with patent infections, needed for prev in 20-25y
  //Nsize is the population size at time t allowing population growth
  int Nsize = (int)(_populationSize * exp(rho*Simulation::simulationTime));
  int survivsSoFar = 0;
  
  // Is the individual in the age range to be pregnant? Set when age reaches appropriate range.
  bool isAtRiskOfFirstPregnancy = false;
  
  int noOfAgeGroupsSharedMem = std::max(Simulation::gMainSummary->getNumOfAgeGroups(),KappaArraySize);
  double* kappaByAge = new double[noOfAgeGroupsSharedMem];
  int* nByAge = new int[noOfAgeGroupsSharedMem];
  for (int i=0; i<noOfAgeGroupsSharedMem; i++) {
    kappaByAge[i] = 0.0;
    nByAge[i] = 0;
  }
  //  Initialise the variable used for calculating infectiousness
  double sumWt_kappa= 0.0;
  double sumWeight  = 0.0;
  
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
    
    //else update the individual
      ++survivsSoFar;
      double ageYears = iter->getAgeInYears(Simulation::simulationTime);
      double availability = iter->infIncidence->_BaselineAvailabilityToMosquitoes * _transmissionModel->getRelativeAvailability(ageYears);
      sumWeight += availability;
      sumWt_kappa += availability*iter->withinHostModel->getProbTransmissionToMosquito();
      
      // kappaByAge and nByAge are used in the screensaver only
      // TODO: This measure of infectiousness isn't directly affected by
      // bed-nets and isn't usable with NC's vector transmission model.
      int ia = iter->ageGroup() - 1;
      kappaByAge[ia] += iter->withinHostModel->getProbTransmissionToMosquito();
      ++nByAge[ia];
      
      /*
      TODO: we need to check the order (ok) that all data
      is in a defined state when used.
      We call outmigrate if this is the last human in the list
      or if their dob differs from the next one (ie last in 5day agegroup)
      */
      //std::cout  << "before om if" <<iter->hData.dob<<" "<<iter->hData.ID<<std::endl;  
      if (iter == last ||	// Last human
          // Copy iter, increment, and compare getDateOfBirth:
          (++HumanIter(iter))->getDateOfBirth() != iter->getDateOfBirth())
        if (outmigrate(*iter, Nsize, survivsSoFar)) {
          iter->destroy();
          iter = _population.erase(iter);
          continue;
        }
    
    //Determine risk from maternal infection from   
    if(ageYears < 25.0) {
      if (ageYears >= 20.0) {
	/* updates the counts of the number of individuals of child bearing age
	and the numbers of these with patent parasitemia */
	
	if(!isAtRiskOfFirstPregnancy) {
	  nCounter = 0;
	  pCounter = 0;
	  isAtRiskOfFirstPregnancy = true;
	}
	nCounter ++;
	if (iter->withinHostModel->parasiteDensityDetectible()){
	  pCounter ++;
	}
      }
      
      if(isAtRiskOfFirstPregnancy && ageYears < 20.0) {
	isAtRiskOfFirstPregnancy = false;	// only call once per time-step
	PathogenesisModel::setRiskFromMaternalInfection(nCounter, pCounter);
      }
    }
    
    ++iter;
  }	// end of per-human updates
  
  // Shared graphics: report infectiousness
  if (Simulation::simulationTime % 6 ==  0) {
    for (int i=0; i < Simulation::gMainSummary->getNumOfAgeGroups(); i++)
      kappaByAge[i] /= nByAge[i];
    SharedGraphics::copyKappa(kappaByAge);
  }
  
  // increase population size to Nsize
  while (survivsSoFar < Nsize) {
    newHuman(Simulation::simulationTime);
    //++nCounter;
    ++survivsSoFar;
  }
  
  // Calculate kappa (total infectiousness)
  // Currently we use the same summed weights as before. Doing them here would
  // be different because of outmigrated individuals and new births.
  _transmissionModel->updateKappa (sumWeight, sumWt_kappa);
  
  delete [] nByAge;
  delete [] kappaByAge;
}

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
    changeHealthSystem (&interv->getChangeHS().get());
    CaseManagementModel::init();	// should re-read all parameters
    
    //TODO: Do we also need to re-init the kappa array?
    _transmissionModel->copyToInitialKappa();
  }
  
  if (interv->getChangeEIR().present()) {
    NonVectorTransmission* nvt = dynamic_cast<NonVectorTransmission*> (_transmissionModel);
    if (nvt != NULL) {
      Global::simulationMode=transientEIRknown;
      nvt->inputEIR(interv->getChangeEIR().get());
    } else {
      throw xml_scenario_error("Warning: changeEIR intervention can only be used with NonVectorTransmission model!");
    }
  }
  
  if (interv->getVaccinate().present()) {
    vaccinatePopulation(interv->getVaccinate().get());
  }
  if (interv->getMDA().present()) {
    massTreatment(interv->getMDA().get());
  }
  if (interv->getIpti().present()) {
    massIPTiTreatment(interv->getIpti().get());
  }
  
  /* TODO
  if (interv->getIRS().present()) {
  } */
}

void Population::massTreatment(const scnXml::Mass& mass){
  double minAge = mass.getMinAge();
  double maxAge = mass.getMaxAge();
  double compliance = mass.getCoverage();
  
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    double ageYears = iter->getAgeInYears(Simulation::simulationTime);
    if ((iter->withinHostModel->getCumulativeInfections() > 0) && (ageYears > minAge) && (ageYears < maxAge)){
      if (W_UNIFORM() < compliance) {
	/* TODO: here we assume a 100% clearance rate for the MDA drug we use.
	This is not consistent with the way we treat according to the Health
	system description. The default clearance rate for MDA should be 100%
	since this simulates what was meant to happen in Garki.  We can change
	this by introducing an optional clearance rate that can be < 100%. */
	iter->clearInfections();
      }
    }
    /* The following line of code is added for calculating expected
    inoculations for the analysis of pre-erythrocytic immunity.
    TODO: inside the above conditional? */
    // Only affects Summary::addToExpectedInfected
    iter->infIncidence->_pinfected = 0.0;
  }
}

void Population::massIPTiTreatment(const scnXml::Mass& mass){
  //Set the last SP Dose given for the eligible humans - is this all we need to do?     
  double minAge = mass.getMinAge();
  double maxAge = mass.getMaxAge();
  double compliance = mass.getCoverage();
  
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter) {
    double ageYears = iter->getAgeInYears(Simulation::simulationTime);
    if ((ageYears > minAge) && (ageYears < maxAge))
      iter->IPTiTreatment(compliance);
  }
}


void Population::vaccinatePopulation(const scnXml::Mass& mass){
  double minAge = mass.getMinAge();
  double maxAge = mass.getMaxAge();
  double compliance = mass.getCoverage();

  for(HumanIter iter=_population.begin(); iter != _population.end(); ++iter){
    double ageYears = iter->getAgeInYears(Simulation::simulationTime);
    if ((ageYears > minAge) && (ageYears < maxAge)) {
      if (W_UNIFORM() < compliance){
        iter->vaccinate();
        Simulation::gMainSummary->reportMassVaccination(iter->ageGroup());
      }
    }
  }
}


short Population::outmigrate(Human& current, int Nsize, int &survivsSoFar){
  // maximum remaining lifespan in timesteps:
  int j=_maxTimestepsPerLife-(Simulation::simulationTime-current.getDateOfBirth());
  
  double targetPop=cumpc[j-1] * Nsize;	//target population in age group of current human
  
  // Actual number of people so far = Survivsofar
  // Number to be removed is the difference, rounded to the nearest integer
  int outmigrs = survivsSoFar - (int)floor(targetPop+0.5);
  // We can't outmigrate more than one person at once:
  if (outmigrs > 1) outmigrs = 1;
  if (outmigrs == 1){
    --survivsSoFar;
    return true;
  }
  return false;
}

// Static method used by estimateRemovalRates
double Population::setDemoParameters (double param1, double param2) {
  rho = get_growthrate() / (100.0*(Global::intervalsPerYear));
  
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
