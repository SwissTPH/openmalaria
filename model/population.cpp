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
#include "caseManagement.h"
#include "intervention.h"
#include "transmissionModel.h"
#include "summary.h"
#include "boincWrapper.h"
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

Population::Population(){

  this->init();

}

Population::Population(TransmissionModel* transmissionModel, int populationSize) 
  : _populationSize(populationSize), _transmissionModel(transmissionModel) {
  this->init();
}

Population::~Population() {
  delete _caseManagement; 
}

void Population::init(){
  _caseManagement = new CaseManagementModel();
  _workUnitIdentifier=get_wu_id();
  _maxTimestepsPerLife=maxLifetimeDays/Global::interval;
  cumpc = (double*)malloc(((_maxTimestepsPerLife))*sizeof(double));

  IDCounter=0;
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
  const AgeGroupPerC::GroupSequence& group = getDemography().getAgeGroup().getGroup();
  if (group.size() < ngroups-1) {
    cerr << "expected " << ngroups-1 << " elements of \"group\" in demography->ageGroup (in scenario.xml)" << endl;
    throw 0;
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
  int predpercX;
  double M1s;
  double M2s;
  double Ms;
  double* predperc = new double[_maxTimestepsPerLife];
  predpercX = _maxTimestepsPerLife;
  cumpc[0] = 0.0;
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    double ageYears = (_maxTimestepsPerLife-j-1) / (1.0*Global::intervalsPerYear);
    M1s=(mu0 * (1.0-exp(-alpha0*ageYears)) / alpha0);
    M2s=(mu1 * (exp(alpha1*ageYears)-1.0) / alpha1);
    Ms=M1s+M2s;
    predperc[j]=exp(-rho*ageYears-Ms);
    if (j < _maxTimestepsPerLife-Global::maxAgeIntervals){
      predperc[j]=0.0;
    }
    cumpc[j]=cumpc[j-1]+predperc[j];
  }
  int N_new=0;
  int cumulativePop=0;
  for (int j=1;j<_maxTimestepsPerLife; j++) {
    int iage=_maxTimestepsPerLife-j-1;
    //Scale using the total cumpc
    cumpc[j]=cumpc[j]/cumpc[_maxTimestepsPerLife-1];
    N_new=(int)floor(cumpc[j]*_populationSize+0.5)-cumulativePop;
    for (int j1=1;j1<=N_new; j1++){
      if (!isCheckpoint){
        newHuman(-iage);
      }
      cumulativePop=cumulativePop+1;
    }
  }
  delete [] predperc;
}

void Population::initialiseInfantArrays(){
  for ( int i=0;i<Global::intervalsPerYear; i++) {
    Global::infantIntervalsAtRisk[i]=0;
    Global::infantDeaths[i]=0;
  }
}

void Population::updateInfantArrays(int agetstep, int doomed){
  ++Global::infantIntervalsAtRisk[agetstep];
  if ((doomed == 4) || (doomed == -6) || (doomed == 6)){
    ++Global::infantDeaths[agetstep];
  }
}

void Population::initialiseHumanList(){
  IDCounter=0;
}

void Population::writeLists (fstream& funit) {
  funit << _transmissionModel->annualEIR << endl;
  funit << _populationSize << endl;
  funit << IDCounter << endl;
  funit << mu0 << endl;
  funit << mu1 << endl;
  funit << alpha0 << endl;
  funit << alpha1 << endl ;
  funit << rho << endl; 
  funit << _annualAverageKappa << endl;
  funit << _sumAnnualKappa << endl;
  for (int i = 0; i < Global::intervalsPerYear; ++i)
    funit << _transmissionModel->kappa[i] << endl;
  funit <<  _workUnitIdentifier << endl;

  //Write human data
  start_cp_timer();
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    funit << *iter;
  }

  //Finished writing lists
  stop_cp_timer();

}

void Population::readLists (fstream& funit) {
  //Start reading a checkpoint
  initialiseHumanList();
  funit >> _transmissionModel->annualEIR;
  funit >> _populationSize;
  funit >> IDCounter;
  funit >> mu0;
  funit >> mu1;
  funit >> alpha0;
  funit >> alpha1;
  funit >> rho;
  funit >> _annualAverageKappa;
  funit >> _sumAnnualKappa;
  for (int i = 0; i < Global::intervalsPerYear; ++i)
    funit >> _transmissionModel->kappa[i];
  funit >> _workUnitIdentifier;

  if ( _workUnitIdentifier !=  get_wu_id()) {
    cout << "cp_ct" << get_wu_id() << _workUnitIdentifier << endl;
    exit(-9);
  }

  //Start reading the human data
  int indCounter = 0;	// Number of individuals read from checkpoint
  while(!(funit.eof()||_populationSize==indCounter)){
      //continue: Fortran cont is probably not C cont
    _population.push_back(Human(funit, _caseManagement, Simulation::simulationTime));
    indCounter++;
  }
  if ((_populationSize !=  get_populationsize()) || (_populationSize !=  indCounter)){
    cout << "exp_p_s" << _populationSize << get_populationsize() << indCounter << endl;
    exit(-7);
  }
}

void Population::newHuman(int dob){
  ++IDCounter;
  _population.push_back(Human(IDCounter, dob,  _caseManagement, Simulation::simulationTime));
}

void Population::update1(){
  //NOTE: this needs to be called somewhere; possibly where its called will affect the results slightly.
  _transmissionModel->advancePeriod (_population, Simulation::simulationTime);
  
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
  int ia = 0;
  
  // Update each human in turn
  //std::cout<<" time " <<t<<std::endl;
  HumanIter last = _population.end();
  --last;
  for (HumanIter iter = _population.begin(); iter != _population.end();){
    int agetstep = Simulation::simulationTime-iter->getDateOfBirth();
    double ageYears = iter->getAgeInYears(Simulation::simulationTime);
    //First eliminate those who die or who are too old
    if (agetstep > Global::maxAgeIntervals) {
      iter->setDoomed(1);
    }
    if (iter->getDoomed() > 0){
      iter->destroy();
      iter=_population.erase(iter);
      continue;
    } else {	//else update the individual
      ++survivsSoFar;
      // UPDATE HUMAN
      iter->update(Simulation::simulationTime,_transmissionModel);
      double availability = iter->getBaselineAvailabilityToMosquitoes() * _transmissionModel->getRelativeAvailability(ageYears);
      sumWeight += availability;
      sumWt_kappa += availability*iter->getProbTransmissionToMosquito();
      
      //  update array for the infant death rates     
      if (agetstep <= Global::intervalsPerYear){
        updateInfantArrays(agetstep-1, iter->getDoomed());
      }
      ia = iter->ageGroup() - 1;
      /*
      TODO: ptransmit should depend on bednet usage
      kappaByAge and nByAge are used in the screensaver only
      */
      kappaByAge[ia] += iter->getProbTransmissionToMosquito();
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
    }
    //Determine risk from maternal infection from   
    if(ageYears < 25.0 && ageYears >= 20.0) {
      updateMaternalMalariaCounters(*iter, isAtRiskOfFirstPregnancy, nCounter, pCounter);
    }
    if(isAtRiskOfFirstPregnancy && ageYears < 20.0) {
      isAtRiskOfFirstPregnancy = false;	// only call once per time-step
      _caseManagement->setRiskFromMaternalInfection(nCounter, pCounter);
    }
    ++iter;
  }
  /*
  end of updating each individual
  now update the population-level measures, and add new births
  */
  if (Simulation::simulationTime % 6 ==  0) {
    for (int i=0; i < Simulation::gMainSummary->getNumOfAgeGroups(); i++) {
      kappaByAge[i] /= nByAge[ia];
    }
    SharedGraphics::copyKappa(kappaByAge);
  }
  int nbirths = Nsize - survivsSoFar;
  for (int i=0;i<nbirths; i++) {
    newHuman(Simulation::simulationTime);
    ++nCounter;
  }
  int tmod = (Simulation::simulationTime - 1) % Global::intervalsPerYear;
  //Prevent NaNs
  if (sumWeight == 0) {
    _transmissionModel->kappa[tmod]=0;
    cout << "sW.eq.0" << endl;
  }
  else {
    _transmissionModel->kappa[tmod] = sumWt_kappa / sumWeight;
  }
  
  //Calculate time-weighted average of kappa
  if (tmod == 0) {
    _sumAnnualKappa=0.0;
  }
  _sumAnnualKappa += _transmissionModel->kappa[tmod] * Global::interval * _transmissionModel->EIR[tmod];
  if (tmod + 1 == Global::intervalsPerYear) {
    if ( _transmissionModel->annualEIR ==  0) {
      _annualAverageKappa=0;
      cout << "aE.eq.0" << endl;
    }
    else {
      _annualAverageKappa=_sumAnnualKappa/_transmissionModel->annualEIR;
    }
  }
  delete [] nByAge;
  delete [] kappaByAge;
}

void Population::newSurvey () {
  for(HumanIter iter=_population.begin(); iter != _population.end(); iter++){
    iter->summarize();
  }
  Simulation::gMainSummary->setNumTransmittingHosts(_transmissionModel->kappa[Global::modIntervalsPerYear(Simulation::simulationTime) - 1]);
  Simulation::gMainSummary->setAnnualAverageKappa(_annualAverageKappa);
  Simulation::gMainSummary->incrementSurveyPeriod();

}

void Population::implementIntervention (int time) {
  const Intervention* interv = getInterventionByTime (time);
  if (interv == NULL)
    return;
  
  if (interv->getChangeHS().present()) {
    changeHealthSystem (&interv->getChangeHS().get());
    delete _caseManagement;
    _caseManagement = new CaseManagementModel();
    for(HumanIter iter=_population.begin(); iter != _population.end(); iter++){
      iter->setCaseManagement(_caseManagement);
    }
    
    //TODO: Do we also need to re-init the kappa array?
    memcpy (_transmissionModel->initialKappa, _transmissionModel->kappa, Global::intervalsPerYear * sizeof(*_transmissionModel->kappa));
  }
  
  if (interv->getChangeEIR().present()) {
    changeEntoData (&interv->getChangeEIR().get());
    Global::simulationMode=transientEIRknown;
    _transmissionModel->inputEIR();
  }
  
  if (interv->getVaccinate().present()) {
    vaccinatePopulation(interv->getVaccinate().get(), time);
  }
  if (interv->getMDA().present()) {
    massTreatment(interv->getMDA().get(), time);
  }
  if (interv->getIpti().present()) {
    massIPTiTreatment(interv->getIpti().get(), time);
  }
  
  /* TODO
  if (interv->getIRS().present()) {
  } */
}

void Population::massTreatment(const Mass& mass, int time){
  double minAge = mass.getMinAge();
  double maxAge = mass.getMaxAge();
  double compliance = mass.getCoverage();
  
    /*
  TODO: here we assume a 100% clearance rate for the MDA drug we use. This is not consistent with
  the way we treat according to the Health system description.
  the default clearance rate for MDA should be 100% since this simulates what was meant to happen
  in Garki.  We can change this by introducing an optional clearance rate that can be < 100%
    */
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    double ageYears = iter->getAgeInYears(Simulation::simulationTime);
    if ((iter->getCumulativeInfections() > 0) && (ageYears > minAge) && (ageYears < maxAge)){
      if (W_UNIFORM() < compliance) {
        iter->clearAllInfections();
      }
    }
        /*
    The following line of code is added for calculating expected inoculations for the analysis
    of pre-erythrocytic immunity
    TODO: inside the above conditional?
        */
    iter->setProbabilityOfInfection(0.0);

  }
}

void Population::massIPTiTreatment(const Mass& mass, int time){
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

void Population::clear() {
  HumanIter iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    iter->destroy();
  }
  _population.clear();
}


void Population::vaccinatePopulation(const Mass& mass, int time){
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
  /*
    
  TODO: I had to extract this from update1 because the surrounding conditional
  would not have worked without dupclicating it. However, I don't understand how it
  works. The change did not break the test case, but we need to go through this with AR.
  fix setting of riskfrommaternalinfection
  */

  int j;  // maximum remaining lifespan in timesteps
  int outmigrs;
  double targetPop;
  bool valoutmigrate=false;
  j=_maxTimestepsPerLife-(Simulation::simulationTime-current.getDateOfBirth());
  targetPop=cumpc[j-1] * Nsize;	//target population in age group of current human
    /*
  Actual number of people so far = Survivsofar
  Number to be removed is the difference, rounded to the nearest integer
    */
  outmigrs = survivsSoFar - (int)floor(targetPop+0.5);
  //std::cout <<" outmigrs "<<outmigrs<<std::endl;
  //Outmigrs should not excced 1, otherwise the calling routine has trouble
  outmigrs=min(outmigrs,1);
  if (outmigrs >= 1){
    --survivsSoFar;
    //std::cout <<" outmigrate "<<current.hData.dob<<"
  //"<<current.hData.ID<<std::endl;
    valoutmigrate=true;
  }
  return valoutmigrate;
}

void Population::updateMaternalMalariaCounters(Human& current, bool &isAtRiskOfFirstPregnancy, int &nCounter, int &pCounter){
  /* updates the counts of the number of individuals of child bearing age
     and the numbers of these with patent parasitemia */
  
  if(!isAtRiskOfFirstPregnancy) {
    nCounter = 0;
    pCounter = 0;
    isAtRiskOfFirstPregnancy = true;
  }
  nCounter ++;
  if (current.getTotalDensity() > Human::detectionlimit){
    pCounter ++;
  }
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
