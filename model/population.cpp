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
#include "openMalaria.h"
#include "inputData.h"
#include "ShmStruct.h"
#include "human.h"
#include "simulation.h"
#include "caseManagement.h"
#include "intervention.h"
#include "transmissionModel.h"
#include "summary.h"
#include <math.h>

using namespace std;

double a0[ngroups];
double a1[ngroups];
double perc[ngroups];
double M[ngroups];
double M1[ngroups];
double M2[ngroups];
double pred[ngroups];
double mu0;
double mu1;
double alpha0;
double alpha1;
double rho;
double *cumpc;
int cumpcX;

int IDCounter;

Population::Population(){

  this->init();

}

Population::Population(TransmissionModel* transmissionModel, int populationSize) 
  : _populationSize(populationSize), _transmissionModel(transmissionModel) {
  this->init();
}

Population::~Population(){

 delete _caseManagement; 
 delete _transmissionModel;  

}

void Population::init(){

  matArray.resize(30);
  _caseManagement = new CaseManagement();
  _workUnitIdentifier=get_wu_id();
  _maxTimestepsPerLife=maxLifetimeDays/interval;
  cumpc = (double*)malloc(((_maxTimestepsPerLife))*sizeof(double));
  cumpcX = _maxTimestepsPerLife;

  IDCounter=0;

}

void Population::estimateRemovalRates () {
   /*  
    
        mu1, alpha1: These are estimated here
        mu0, alpha0: These are fixed alpha0=4.0, mu0 calculated from the other parameters
        rho - population growth rate (input)
    
    */
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
    for (int i=1;i<ngroups; i++) {
      a1[i]=get_demo_upperbound(i-1);
      perc[i]=get_popperc(i-1);
      sumperc=sumperc+perc[i];
    }
    for (int i=2;i<ngroups; i++) {
      a0[i]=get_demo_upperbound(i-2);
    }
    //Add age group for first month of life
    a0[0]=0.0;
    a1[0]=1.0/12.0;
    a0[1]=a1[1 - 1];
    perc[0]=0.0;
    for (int i=0;i<ngroups; i++){
      perc[i]=perc[i]/sumperc*100;
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
  double yage;
  double *predperc;
  int predpercX;
  double M1s;
  double M2s;
  double Ms;
  int iage;
  predperc = new double[_maxTimestepsPerLife];
  predpercX = _maxTimestepsPerLife;
  cumpc[0]=0.0;
  for (int j=2;j<=_maxTimestepsPerLife; j++) {
    iage=_maxTimestepsPerLife-j;
    yage=(_maxTimestepsPerLife-j)/(1.0*intervalsPerYear);
    M1s=(mu0*(1-exp(-alpha0*yage))/alpha0);
    M2s=(mu1*(exp(alpha1*yage)-1)/alpha1);
    Ms=M1s+M2s;
    predperc[j - 1]=exp(-rho*yage-Ms);
    if (j <= _maxTimestepsPerLife-maxAgeIntervals){
      predperc[j-1]=0.0;
    }
    cumpc[j-1]=cumpc[j-1-1]+predperc[j-1];
  }
  int N_new=0;
  int cumulativePop=0;
  for (int j=2;j<=_maxTimestepsPerLife; j++) {
    iage=_maxTimestepsPerLife-j;
    //Scale using the total cumpc
    cumpc[j-1]=cumpc[j-1]/cumpc[_maxTimestepsPerLife-1];
    N_new=(int)floor(cumpc[j-1]*_populationSize+0.5)-cumulativePop;
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
    /*
    
    init infant arrays used for infant mortality rates 
    
    */
  for ( int i=0;i<intervalsPerYear; i++) {
    infantIntervalsAtRisk[i]=0;
    infantDeaths[i]=0;
  }
}

void Population::updateInfantArrays(int agetstep, int doomed){
  infantIntervalsAtRisk[agetstep-1]=infantIntervalsAtRisk[agetstep-1]+1;
  if (( doomed ==  4) || ( doomed ==  -6) || ( doomed ==  6)){
    infantDeaths[agetstep - 1]=infantDeaths[agetstep-1]+1;
  }
}

void Population::initialiseHumanList(){
  IDCounter=0;
}

void Population::writeLists (fstream& funit) {
   
  /*
    Start write checkpoint
    Write parameters
  */
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
  funit << *_transmissionModel->kappa << endl ;
  funit <<  _workUnitIdentifier << endl;

  //Write human data
  start_cp_timer();
  std::list<Human*>::iterator iter;
  for(iter=_population.begin(); iter != _population.end(); ++iter){
    funit << **iter;
  }

  //Finished writing lists
  stop_cp_timer();

}

void Population::readLists (fstream& funit) {
   
    //int i;
    //Counter for number of individuals stored in the checkpoint
    int indCounter;
    //Start reading a checkpoint
    initialiseHumanList();
    indCounter=0;
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
    funit >> *_transmissionModel->kappa; 
    funit >> _workUnitIdentifier;

    if ( _workUnitIdentifier !=  get_wu_id()) {
      cout << "cp_ct" << get_wu_id() << _workUnitIdentifier << endl;
      exit(-9);
    }

    //Start reading the human data
    while(!(funit.eof()||_populationSize==indCounter)){
      //continue: Fortran cont is probably not C cont
      _population.push_back(new Human(funit, _caseManagement, Simulation::simulationTime));
      indCounter++;
    }
    if ((_populationSize !=  get_populationsize()) || (_populationSize !=  indCounter)){
      cout << "exp_p_s" << _populationSize << get_populationsize() << indCounter << endl;
      exit(-7);
    }
}

void Population::newHuman(int dob){
    IDCounter=IDCounter+1;
    // TO DO Check memory leaks
    _population.push_back(new Human(IDCounter, dob,  _caseManagement, Simulation::simulationTime));
    //delete human;

}

void Population::update1(){
  /*
    
    Updates all individuals in the list for one time-step
    
    Also updates the population-level measures such as infectiousness, and the age-distribution by c
    outmigrating or creating new births if necessary
    
  */

  int nbirths;
  int tmod;
  double availability;
  //NCounter is the number of indivs per demogr age group
  int nCounter;
  //PCounter is the number with patent infections, needed for prev in 20-25y
  int pCounter;
  int ia=0;
  //Nsize is the population size at time t allowing population growth 
  int Nsize;
  short currentIsRemoved;
  int survivsSoFar;
  int i;
  int k;
  int agetstep;
  //age in years, used to store age of individual who are removed unit the end of the loop
  double yage;
  double sumWt_kappa;
  double sumWeight;
  double sumAvailability;
  double sumAvailabilityITN;
  double adultEIR;
  int noOfAgeGroupsSharedMem = std::max(Simulation::gMainSummary->getNumOfAgeGroups(),KappaArraySize);
  double* kappaByAge = new double[noOfAgeGroupsSharedMem];
  int* nByAge = new int[noOfAgeGroupsSharedMem];
  k=ngroups-1;
  survivsSoFar=0;
  nCounter=0;
  pCounter=0;
  Nsize=(int)(_populationSize*exp((rho)*Simulation::simulationTime));
  adultEIR=_transmissionModel->calculateEIR(Simulation::simulationTime);
  for ( i=0;i<noOfAgeGroupsSharedMem; i++) {
    kappaByAge[i]=0.0;
    nByAge[i]=0;
  }
  //  Initialise the variable used for calculating infectiousness
  sumWt_kappa=0.0;
  sumWeight=0.0;
  sumAvailability=0.0;
  sumAvailabilityITN=0.0;
  // Update each human in turn
  //std::cout<<" time " <<t<<std::endl;
  std::list<Human*>::iterator end = _population.end(), last = _population.end();
  --last;
  std::list<Human*>::iterator iter= _population.begin();
  while( iter != end ){
        
    currentIsRemoved=0/* bool */;
    agetstep=Simulation::simulationTime-(*iter)->getDateOfBirth();
    yage=(*iter)->getAgeInYears(Simulation::simulationTime);
    //First eliminate those who die or who are too old
    if ( agetstep >  maxAgeIntervals) {
      (*iter)->setDoomed(1);
    }
    if ((*iter)->getDoomed() > 0){
      delete *iter;
      iter=_population.erase(iter);
      currentIsRemoved=1/* bool */;
      //else update the individual 
    }
    else {
      survivsSoFar=survivsSoFar+1;
      nCounter=nCounter+1;
      if ((*iter)->getTotalDensity() > detectionlimit){
        pCounter=pCounter+1;
      }       
      // UPDATE HUMAN
      (*iter)->update(Simulation::simulationTime,_transmissionModel);
      availability=(*iter)->getBaselineAvailabilityToMosquitoes()*_transmissionModel->getRelativeAvailability(yage);
      sumAvailability=sumAvailability+availability;
      //TODO: check ITN code
      if ( (*iter)->hasInsecticideTreatedNet()) {
        sumAvailabilityITN=sumAvailabilityITN+availability;
        sumWeight=sumWeight+availability*Pu1/Pu0;
        sumWt_kappa=sumWt_kappa+availability*(*iter)->getProbTransmissionToMosquito()*Pu1/Pu0;
      }
      else {
        sumWeight=sumWeight+availability;
        sumWt_kappa=sumWt_kappa+availability*(*iter)->getProbTransmissionToMosquito();
      }
      //  update array for the infant death rates     
      if (agetstep <= intervalsPerYear){
        updateInfantArrays(agetstep, (*iter)->getDoomed());
      }
      ia=(*iter)->ageGroup();
      /*
        TODO: ptransmit should depend on bednet usage
        kappaByAge and nByAge are used in the screensaver only
      */
      kappaByAge[ia-1]=kappaByAge[ia-1]+(*iter)->getProbTransmissionToMosquito();
      nByAge[ia-1]=nByAge[ia-1]+1;
        
      /*
        TODO: we need to check the order (ok) that all data
        is in a defined state when used.
        We call outmigrate if this is the last human in the list
        or if their dob differs from the next one (ie last in 5day agegroup)
      */
      //std::cout  << "before om if" <<iter->hData.dob<<" "<<iter->hData.ID<<std::endl;  
      if (iter==last){
        currentIsRemoved=outmigrate(&(**iter), Nsize,&( survivsSoFar),&( nCounter),&( pCounter),&( k), yage);
        if ( currentIsRemoved >  0) {
          delete *iter;
          iter=_population.erase(iter);
          //std::cout << " new list pointer "<<iter->hData.dob<<" "<<iter->hData.ID<<std::endl;
        }

      }
      else{

        //TODO: establish outmigrate requirements, and clean up the following.
        ++iter;
        int dob=(*iter)->getDateOfBirth();
        --iter;
        if ((*iter)->getDateOfBirth() != dob){
          currentIsRemoved=outmigrate((*iter), Nsize,&( survivsSoFar),&( nCounter),&( pCounter),&( k), yage);
          if ( currentIsRemoved >  0) {
            delete *iter;
            iter=_population.erase(iter);
            //std::cout << " new list pointer "<<iter->hData.dob<<" "<<iter->hData.ID<<std::endl;
          }
          
        }
      }
    }
    if ( currentIsRemoved ==  0) {
      ++iter;
    }
  }
  /*
    end of updating each individual
    now update the population-level measures, and add new births
  */
  if ( mymodf(Simulation::simulationTime, 6) ==  0) {
    for ( i=1;i<=Simulation::gMainSummary->getNumOfAgeGroups(); i++) {
      kappaByAge[i - 1]=kappaByAge[i - 1]/nByAge[ia - 1];
    }
    add_kappa(kappaByAge);
  }
  nbirths=Nsize-survivsSoFar;
  for ( i=1;i<=nbirths; i++) {
    newHuman(Simulation::simulationTime);
    nCounter=nCounter+1;
  }
  tmod=modIntervalsPerYear(Simulation::simulationTime);
  //Prevent NaNs
  if ( sumWeight ==  0) {
    _transmissionModel->kappa[tmod - 1]=0;
    cout << "sW.eq.0" << endl;
  }
  else {
    _transmissionModel->kappa[tmod - 1]=sumWt_kappa/sumWeight;
  }
  z=sumAvailabilityITN/sumAvailability;
  //Calculate time-weighted average of kappa
  if ( tmod ==  1) {
    _sumAnnualKappa=0.0;
  }
  _sumAnnualKappa=_sumAnnualKappa+_transmissionModel->kappa[tmod - 1]*interval*_transmissionModel->EIR[tmod - 1];
  if ( tmod ==  intervalsPerYear) {
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
  
  std::list<Human*>::iterator iter;
  for(iter=_population.begin(); iter != _population.end(); iter++){
    (*iter)->summarize();
  }
  Simulation::gMainSummary->setNumTransmittingHosts(_transmissionModel->kappa[modIntervalsPerYear(Simulation::simulationTime) - 1]);
  Simulation::gMainSummary->setAnnualAverageKappa(_annualAverageKappa);
  Simulation::gMainSummary->incrementSurveyPeriod();

}

void Population::implementIntervention (int time) {
    //DOCU
    int currentIntervention;
    int i;
    currentIntervention=get_intervention(time);
    if ( currentIntervention ==  no_intervention) {
        return ;
    }
    if ( isOptionIncluded(currentIntervention, change_hs_intervention)) {
      delete _caseManagement;
      _caseManagement = new CaseManagement();
      std::list<Human*>::iterator iter;
      for(iter=_population.begin(); iter != _population.end(); iter++){
        (*iter)->setCaseManagement(_caseManagement);
      }

      //TODO: Do we also need to re-init the kappa array?
      for ( i=1;i<=intervalsPerYear; i++) {
        _transmissionModel->initialKappa[i - 1]=_transmissionModel->kappa[i - 1];
      }
    }
    if ( isOptionIncluded(currentIntervention, change_eir_intervention)) {
        simulationMode=transientEIRknown;
        _transmissionModel->inputEIR();
    }
    if ( isOptionIncluded(currentIntervention, vaccine_intervention)) {
        vaccinatePopulation(time);
    }
    if ( isOptionIncluded(currentIntervention, mda_intervention)) {
        massTreatment(time);
    }
    if ( isOptionIncluded(currentIntervention, ipti_intervention)) {
        massIPTiTreatment(time);
    }
    /*
    TODO: This is meant to be IRS.
    Probably need this to model vector control in scenarios without intervention field ento data.
    if (isOptionIncluded(currentIntervention,irs_intervention)) then
      call vectorControl(time)
    endif
    */
}

void Population::massTreatment(int time){
   
    double maxAge;
    double minAge;
    double compliance;
    double ageYears;
    minAge=get_minage_mda(time);
    maxAge=get_maxage_mda(time);
    /*
    TODO: here we assume a 100% clearance rate for the MDA drug we use. This is not consistent with
    the way we treat according to the Health system description.
    the default clearance rate for MDA should be 100% since this simulates what was meant to happen
    in Garki.  We can change this by introducing an optional clearance rate that can be < 100%
    */
    compliance=get_coverage_mda(time);
    std::list<Human*>::iterator iter;
    for(iter=_population.begin(); iter != _population.end(); ++iter){
      ageYears=(*iter)->getAgeInYears(Simulation::simulationTime);
      if (( (*iter)->getCumulativeInfections() > 0) && ( ageYears > minAge) && ( ageYears < maxAge)){
          if ( (W_UNIFORM()) < compliance) {
            (*iter)->clearAllInfections();
          }
        }
        /*
          The following line of code is added for calculating expected inoculations for the analysis
          of pre-erythrocytic immunity
          TODO: inside the above conditional?
        */
      (*iter)->setProbabilityOfInfection(0.0);

    }
}

void Population::massIPTiTreatment(int time){
    double maxAge;
    double minAge;
    double compliance;
    double ageYears;
    minAge=get_minage_ipti(time);
    maxAge=get_maxage_ipti(time);
    //Set the last SP Dose given for the eligible humans - is this all we need to do?     
    compliance=get_coverage_ipti(time);
    std::list<Human*>::iterator iter;
    for(iter=_population.begin(); iter != _population.end(); ++iter){
      ageYears=(*iter)->getAgeInYears(Simulation::simulationTime);
      if (((*iter)->getCumulativeInfections() > 0) && ( ageYears > minAge) && ( ageYears < maxAge)){
        if ((W_UNIFORM()) < compliance){
          (*iter)->setLastIPTIorPlacebo(Simulation::simulationTime);
          /*
 *            iptiEffect denotes treatment or placebo group
 *                       and also the treatment given when sick (trial-dependent)
 *                                 */
          if ( iptiEffect >= 10){
            (*iter)->setSPDose(Simulation::simulationTime);
            Simulation::gMainSummary->reportIPTDose((*iter)->ageGroup());
          }
        }
      }
  }
}

void Population::clear() {
   std::list<Human*>::iterator iter;
   for(iter=_population.begin(); iter != _population.end(); ++iter)
     delete *iter;
  _population.clear();
}


void Population::vaccinatePopulation(int time){
 
    double minAge;
    double maxAge;
    double coverage;
    double ageYears;
    minAge=get_minage_vaccine(time);
    maxAge=get_maxage_vaccine(time);
    coverage=get_coverage_mass_vaccine(time);
    std::list<Human*>::iterator iter;
    for(iter=_population.begin(); iter != _population.end(); ++iter){
      ageYears=(*iter)->getAgeInYears(Simulation::simulationTime);
      if (( ageYears > minAge) && ( ageYears < maxAge)) {
        if ((W_UNIFORM()) < coverage){
          (*iter)->vaccinate();
          Simulation::gMainSummary->reportMassVaccination((*iter)->ageGroup());
        }
      }
    }
}


short Population::outmigrate(Human *current, int Nsize, int *survivsSoFar, int *nCounter, int *pCounter, int *k, double yage){
  /*
    
    TODO: I had to extract this from update1 because the surrounding conditional
    would not have worked without dupclicating it. However, I don't understand how it
    works. The change did not break the test case, but we need to go through this with AR.
    fix setting of riskfrommaternalinfection
    */

    //k is the number of demogr age groups  
    //Goodman estimated for neonatal mortality due to malaria in pregnancy
    const double gEst= 0.011;
    //Critical value of Prev20-25 for neonatal mortality
    const double critPrev2025= 0.25;
    //Critical value for estimating prevalence in primigravidae
    const double critPrevPrim= 0.19;
    //Proportion of births with primigravid mothers
    const double pBirthPrim= 0.3;
    int h;
    int i;
    int j;
    int outmigrs;
    double targetPop;
    double lbound;
    double prev2025;
    double maxprev;
    //double awidth;
    double prevpg;
    short valoutmigrate;
    valoutmigrate=0;/* bool */
    j=_maxTimestepsPerLife-(Simulation::simulationTime-current->getDateOfBirth());
    targetPop=cumpc[j - 1]*Nsize;
    /*
         Actual number of people so far = Survivsofar
         Number to be removed is the difference, rounded to the nearest integer
    */
    outmigrs=*survivsSoFar-(int)floor(targetPop+0.5);
  //std::cout <<" outmigrs "<<outmigrs<<std::endl;
  //Outmigrs should not excced 1, otherwise the calling routine has trouble
  outmigrs=min(outmigrs,1);
    for (i=1;i<=outmigrs; i++){
      *survivsSoFar =*survivsSoFar-1;
      nCounter =nCounter-1;
      if (current->getTotalDensity() > detectionlimit){
        pCounter =pCounter-1;
      }
  //std::cout <<" outmigrate "<<current->hData.dob<<"
  //"<<current->hData.ID<<std::endl;
      valoutmigrate=1/* bool */;
    }
    /*
      Counts up the number of individuals in each demography age group
      tstep is time from the end of warm up period, t is start of simulation
    */
    if (*k ==1){
      lbound=(get_demo_lowerbound());
    }
    else{
      lbound=get_demo_upperbound(*k-2);
    }
    //  default value for prev2025, for short simulations 
    prev2025=0.25;
    if (yage <  lbound){

        if ( lbound ==  20) {
            prev2025= 1.0 * *pCounter / *nCounter;
            maxprev=prev2025;
            for ( h=2;h<=30; h++) {
                matArray[h-1 - 1]=matArray[h - 1];
                if ( matArray[h-1 - 1] >  maxprev) {
                    maxprev=matArray[h-1 - 1];
                }
            }
            matArray[30 - 1]=prev2025;
            prevpg=1-(1/(1+(maxprev/critPrevPrim)));
            _caseManagement->setRiskFromMaternalInfection(gEst*pBirthPrim*(1-exp(-prevpg/critPrev2025)));
        }
        *k = *k - 1;
        pCounter =0;
        nCounter =0;
    }
    return valoutmigrate;
}


double Population::setDemoParameters (double param1, double param2) {

  return 0;
}
