/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School
 Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
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

#include "summary.h"
#include "inputData.h"

#include "Clinical/event.h"
#include "BoincWrapper.h"

Summary::Summary() : _surveyPeriod(0) {
  //Read all summary parameters from the input file and allocate datastructures.
  const scnXml::AgeGroup::GroupSequence& groups = getMonitoring().getAgeGroup().getGroup();
  _numOfAgeGroups = groups.size()+1;
  _upperbound.resize(_numOfAgeGroups);
  _lowerbound=get_lowerbound();
  /*
    note that the last age group includes individuals who are
    either younger than Lowerbound or older than the last Upperbound
  */

  for (int i=0;i<_numOfAgeGroups-1; i++) {
    _upperbound[i] = groups[i].getUpperbound();
  }
  _upperbound[_numOfAgeGroups-1]=9999;
  _nonMalariaMortality=getParameter(Params::NON_MALARIA_INFANT_MORTALITY);

}

void Summary::initialiseSummaries () {
/*
    
    Allocates memory for event and parasitological summaries. We can wait 
    with the call to this until the beginning of the main simulation.
   
  */
 int numberOfSurveys;
  int i;
  _surveyPeriod=1;

  numberOfSurveys=get_number_of_surveys();
  _surveysTimeIntervals.resize(numberOfSurveys+1);

  for (i=0;i<numberOfSurveys; i++) {
    _surveysTimeIntervals[i]=get_time_of_survey(i);
  }

  _summaryDimensionX=numberOfSurveys+1;
  _summaryDimensionY=_numOfAgeGroups;
  _summaryOption=get_summary_option();
  _assimilatorMode=get_assim_mode();

  /* Vectors initialization */ 
  _numHosts.resize(numberOfSurveys+1);
  _numInfectedHosts.resize(numberOfSurveys+1);
  _numExpectedInfected.resize(numberOfSurveys+1);
  _numPatentHosts.resize(numberOfSurveys+1);
  _sumX.resize(numberOfSurveys+1);
  _sumLogDensity.resize(numberOfSurveys+1);
  _totalInfections.resize(numberOfSurveys+1);
  _numTransmittingHosts.resize(numberOfSurveys+1);
  _totalPatentInfections.resize(numberOfSurveys+1);
  _contributionImmunity.resize(numberOfSurveys+1);
  _pyrogenicThreshold.resize(numberOfSurveys+1);
  _numTreatments1.resize(numberOfSurveys+1);
  _numTreatments2.resize(numberOfSurveys+1);
  _numTreatments3.resize(numberOfSurveys+1);
  _numNonMalariaFever.resize(numberOfSurveys+1);
  _numUncomplicatedEpisodes.resize(numberOfSurveys+1);
  _numSevereEpisodes.resize(numberOfSurveys+1);
  _numDirectDeaths.resize(numberOfSurveys+1);
  _numIndirectDeaths.resize(numberOfSurveys+1);
  _numSequelae.resize(numberOfSurveys+1);
  _numHospitalDeaths.resize(numberOfSurveys+1);
  _numHospitalRecoveries.resize(numberOfSurveys+1);
  _numHospitalSequelae.resize(numberOfSurveys+1);
  _numEPIVaccines.resize(numberOfSurveys+1);
  _numMassVaccines.resize(numberOfSurveys+1);
  _numIPTDoses.resize(numberOfSurveys+1);
  _annualAverageKappa.resize(numberOfSurveys+1);

  for (int i=0;i<numberOfSurveys+1;i++){

    _numHosts[i].resize(_numOfAgeGroups);
    _numInfectedHosts[i].resize(_numOfAgeGroups);
    _numExpectedInfected[i].resize(_numOfAgeGroups);
    _numPatentHosts[i].resize(_numOfAgeGroups);
    _sumX[i].resize(_numOfAgeGroups);
    _sumLogDensity[i].resize(_numOfAgeGroups);
    _totalInfections[i].resize(_numOfAgeGroups);
    _totalPatentInfections[i].resize(_numOfAgeGroups);
    _contributionImmunity[i].resize(_numOfAgeGroups);
    _pyrogenicThreshold[i].resize(_numOfAgeGroups);
    _numTreatments1[i].resize(_numOfAgeGroups);
    _numTreatments2[i].resize(_numOfAgeGroups);
    _numTreatments3[i].resize(_numOfAgeGroups);
    _numNonMalariaFever[i].resize(_numOfAgeGroups);
    _numUncomplicatedEpisodes[i].resize(_numOfAgeGroups);
    _numSevereEpisodes[i].resize(_numOfAgeGroups);
    _numDirectDeaths[i].resize(_numOfAgeGroups);
    _numIndirectDeaths[i].resize(_numOfAgeGroups);
    _numSequelae[i].resize(_numOfAgeGroups);
    _numHospitalDeaths[i].resize(_numOfAgeGroups);
    _numHospitalRecoveries[i].resize(_numOfAgeGroups);
    _numHospitalSequelae[i].resize(_numOfAgeGroups);
    _numEPIVaccines[i].resize(_numOfAgeGroups);
    _numMassVaccines[i].resize(_numOfAgeGroups);
    _numIPTDoses[i].resize(_numOfAgeGroups);   

  }
 


}

Summary::~Summary(){}

void Summary::clearSummaryParameters () {
  _upperbound.clear();
}


void Summary::report(Event& event){
  int reportAgeGroup = event.getAgeGroup() - 1;
  int surveyPeriod = event.getSurveyPeriod() - 1;

  //No reporting during warmup
  if (surveyPeriod < 0)
    return;

  switch (event.getDiagnosis()) {
  case Diagnosis::NON_MALARIA_FEVER:
    _numNonMalariaFever[surveyPeriod][reportAgeGroup]++;
    break;
  case Diagnosis::UNCOMPLICATED_MALARIA:
    _numUncomplicatedEpisodes[surveyPeriod][reportAgeGroup]++;
    break;
  case Diagnosis::SEVERE_MALARIA:
    _numSevereEpisodes[surveyPeriod][reportAgeGroup]++;
    break;
  case Diagnosis::INDIRECT_MALARIA_DEATH:
    _numIndirectDeaths[surveyPeriod][reportAgeGroup]++;
    break;
  default:
    //Diagnosis not conclusive
    cout << "diag nc or non-malaria fever" << lineEnd;
    break;
  }
  switch (event.getOutcome()) {
  case Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_NON_TREATED:
    //do nothing
    break;
  case Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_OUTPATIENTS:
    //do nothing
    break;
  case Outcome::NO_CHANGE_IN_PARASITOLOGICAL_STATUS_INPATIENTS:
    //do nothing
    break;
  case Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_NON_TREATED:
    //do nothing
    break;
  case Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_OUTPATIENTS:
    //do nothing
    break;
  case Outcome::PARASITES_ARE_CLEARED_PATIENT_RECOVERS_INPATIENTS:
    _numHospitalRecoveries[surveyPeriod][reportAgeGroup]++;
    //TODO: we curr do not distinquish between treated vs untreated seqs.
    break;
  case Outcome::PARASITES_ARE_CLEARED_PATIENT_HAS_SEQUELAE_NON_TREATED:
    _numSequelae[surveyPeriod][reportAgeGroup]++;
    break;
  case Outcome::PARASITES_ARE_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS:
    _numSequelae[surveyPeriod][reportAgeGroup]++;
    _numHospitalSequelae[surveyPeriod][reportAgeGroup]++;
    break;
  case Outcome::PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_NON_TREATED:
    _numSequelae[surveyPeriod][reportAgeGroup]++;
    break;
  case Outcome::PARASITES_NOT_CLEARED_PATIENT_HAS_SEQUELAE_INPATIENTS:
    _numSequelae[surveyPeriod][reportAgeGroup]++;
    _numHospitalSequelae[surveyPeriod][reportAgeGroup]++;
    break;
  case Outcome::PATIENT_DIES_NON_TREATED:
    _numDirectDeaths[surveyPeriod][reportAgeGroup]++;
    break;
  case Outcome::PATIENT_DIES_INPATIENTS:
    _numDirectDeaths[surveyPeriod][reportAgeGroup]++;
    _numHospitalDeaths[surveyPeriod][reportAgeGroup]++;
    break;
  case Outcome::INDIRECT_DEATH:
    //do nothing
    break;
  default:
    //Outcome not conclusive
    cout << "oc nc" << event.getOutcome() << lineEnd;
    break;
  }
}

void Summary::report(Pathogenesis::State state, int ageGroup,int surveyPeriod) {
  //No reporting during warmup
  if (surveyPeriod < 1)
    return;
  
  //NOTE: this happens slightly differently to old reporting; we still report
  // sickness when the patient has indirect death.
  if (state & Pathogenesis::MALARIA) {
    if (state & Pathogenesis::COMPLICATED)
      _numSevereEpisodes[surveyPeriod-1][ageGroup-1]++;
    else // UC or UC2
      _numUncomplicatedEpisodes[surveyPeriod-1][ageGroup-1]++;
  } else if (state & Pathogenesis::SICK) {
    _numNonMalariaFever[surveyPeriod-1][ageGroup-1]++;
  }	// also possibility of nothing, but not reported in this case
  
  if (state & Pathogenesis::DIRECT_DEATH)
    _numDirectDeaths[surveyPeriod-1][ageGroup-1]++;
  else if (state & Pathogenesis::INDIRECT_MORTALITY)
    _numIndirectDeaths[surveyPeriod-1][ageGroup-1]++;
  else if (state & Pathogenesis::SEQUELAE)
    _numSequelae[surveyPeriod-1][ageGroup-1]++;
  
  //TODO: we don't know whether patient was in hospital or not
  //_numHospitalRecoveries[surveyPeriod-1][ageGroup-1]++;
  //_numHospitalSequelae[surveyPeriod-1][ageGroup-1]++;
  //_numHospitalDeaths[surveyPeriod-1][ageGroup-1]++;
}

void Summary::reportTreatment (int ageGroup, int regimen) {
    
    //No reporting during warmup
    if ( _surveyPeriod <  1) 
      return;

    switch (regimen) {
    case 1:
      _numTreatments1[_surveyPeriod-1][ageGroup - 1]++;
      break;
    case 2:
      _numTreatments2[_surveyPeriod-1][ageGroup - 1]++;
      break;
    case 3:
      _numTreatments3[_surveyPeriod-1][ageGroup - 1]++;
      break;
    default:
      //Regimen not conclusive
      cout << "reg nc" << regimen << lineEnd;
      break;
    }
}

void Summary::reportEPIVaccination (int ageGroup) {
   
    //No reporting during warmup
    if ( _surveyPeriod <  1) 
      return ;
    
    _numEPIVaccines[_surveyPeriod-1][ageGroup - 1]++;
}

void Summary::reportMassVaccination (int ageGroup) {
    
  //No reporting during warmup
  if ( _surveyPeriod <  1) 
    return;
  
  _numMassVaccines[_surveyPeriod-1][ageGroup - 1]++;

}

void Summary::reportIPTDose (int ageGroup) {
  
    //No reporting during warmup
    if ( _surveyPeriod <  1) 
      return;

    _numIPTDoses[_surveyPeriod-1][ageGroup - 1]++;
}

int Summary::ageGroup(double age){
    if (age < _lowerbound)
      return _numOfAgeGroups;
    
    int valageGroup=1;
    while( age > _upperbound[valageGroup - 1]) {
      valageGroup++;
    }
    return valageGroup;
}

void Summary::printHosts(ostream& out) {

  for(int j=0; j< (int) _numHosts[0].size()-1; j++){

    if (_assimilatorMode ==  1)
      cout << "\t" << _numHosts[0][j] << lineEnd;
    else
      cout << "\t" << 1 << "\t" << j+1 << "\t" 
         << nHost << "\t" <<  _numHosts[0][j] <<  lineEnd;
    
  }

}

void Summary::writeSummaryArrays () {
  string output_filename =
      BoincWrapper::resolveFile("output.txt");
  ifstream test (output_filename.c_str());
  if (test.is_open())
    throw runtime_error ("File output.txt exists!");
  
  ofstream outputFile;
  outputFile.open(output_filename.c_str());
  
  outputFile.width (0);
  // For additional control:
//   outputFile.precision (6);
//   outputFile << scientific;

  if ( isOptionIncluded(_summaryOption, nHost)) {
    writeArray(outputFile, nHost, _assimilatorMode,  _numHosts);
  }
  if ( isOptionIncluded(_summaryOption, nInfect)) {
    writeArray(outputFile, nInfect, _assimilatorMode, _numInfectedHosts);
  }
  if ( isOptionIncluded(_summaryOption, nExpectd)) {
    writeArray(outputFile, nExpectd, _assimilatorMode, _numExpectedInfected);
  }
  if ( isOptionIncluded(_summaryOption, nPatent)) {
    writeArray(outputFile, nPatent, _assimilatorMode, _numPatentHosts);
  }
  if ( isOptionIncluded(_summaryOption, sumX)) {
    writeArray(outputFile, sumX, _assimilatorMode, _sumX);
  }
  if ( isOptionIncluded(_summaryOption, sumlogDens)) {
    writeArray(outputFile, sumlogDens, _assimilatorMode, _sumLogDensity);
  }
  if ( isOptionIncluded(_summaryOption, totalInfs)) {
    writeArray(outputFile, totalInfs, _assimilatorMode, _totalInfections);
  }
  if ( isOptionIncluded(_summaryOption, nTransmit)) {
    writeArray(outputFile, nTransmit, _assimilatorMode, _numTransmittingHosts);
  }
  if ( isOptionIncluded(_summaryOption, totalPatentInf)) {
    writeArray(outputFile, totalPatentInf, 
               _assimilatorMode, _totalPatentInfections);
  }
  if ( isOptionIncluded(_summaryOption, contrib)) {
    writeArray(outputFile, contrib, _assimilatorMode, _contributionImmunity);
  }
  if ( isOptionIncluded(_summaryOption, pyrogenThrs)) {
    writeArray(outputFile, pyrogenThrs, _assimilatorMode, _pyrogenicThreshold);
  }
  if ( isOptionIncluded(_summaryOption, nTreatments1)) {
    writeArray(outputFile, nTreatments1, _assimilatorMode, _numTreatments1);
  }
  if ( isOptionIncluded(_summaryOption, nTreatments2)) {
    writeArray(outputFile, nTreatments2, _assimilatorMode, _numTreatments2);
  }
  if ( isOptionIncluded(_summaryOption, nTreatments3)) {
    writeArray(outputFile, nTreatments3, _assimilatorMode, _numTreatments3);
  }
  if ( isOptionIncluded(_summaryOption, nUncomp)) {
    writeArray(outputFile, nUncomp, 
               _assimilatorMode, _numUncomplicatedEpisodes);
  }
  if ( isOptionIncluded(_summaryOption, nSevere)) {
    writeArray(outputFile, nSevere, _assimilatorMode, _numSevereEpisodes);
  }
  if ( isOptionIncluded(_summaryOption, nSeq)) {
    writeArray(outputFile, nSeq, _assimilatorMode, _numSequelae);
  }
  if ( isOptionIncluded(_summaryOption, nHospitalDeaths)) {
    writeArray(outputFile, nHospitalDeaths, 
               _assimilatorMode, _numHospitalDeaths);
  }
  if ( isOptionIncluded(_summaryOption, nIndDeaths)) {
    writeArray(outputFile, nIndDeaths, _assimilatorMode, _numIndirectDeaths);
  }
  if ( isOptionIncluded(_summaryOption, nDirDeaths)) {
    writeArray(outputFile, nDirDeaths, _assimilatorMode, _numDirectDeaths);
  }
  if ( isOptionIncluded(_summaryOption, nEPIVaccines)) {
    writeArray(outputFile, nEPIVaccines, _assimilatorMode, _numEPIVaccines);
  }
  if ( isOptionIncluded(_summaryOption, nMassVaccines)) {
    writeArray(outputFile, nMassVaccines, _assimilatorMode, _numMassVaccines);
  }
  if ( isOptionIncluded(_summaryOption, nHospitalRecovs)) {
    writeArray(outputFile, nHospitalRecovs, 
               _assimilatorMode, _numHospitalRecoveries);
  }
  if ( isOptionIncluded(_summaryOption, nHospitalSeqs)) {
    writeArray(outputFile, nHospitalSeqs, 
               _assimilatorMode, _numHospitalSequelae);
  }
  if ( isOptionIncluded(_summaryOption, nIPTDoses)) {
    writeArray(outputFile, nIPTDoses, _assimilatorMode, _numIPTDoses);
  }
  if ( isOptionIncluded(_summaryOption, annAvgK)) {
    writeArray(outputFile, annAvgK, _assimilatorMode, _annualAverageKappa);
  }
  if ( isOptionIncluded(_summaryOption, nNMFever)) {
    writeArray(outputFile, nNMFever, _assimilatorMode, _numNonMalariaFever);
  }
  //Infant mortality rate is a single number, therefore treated separately
  if ( isOptionIncluded(_summaryOption, imr_summary)) {
    outputFile << "\t" << 1 << "\t" << 1 << "\t" << imr_summary 
               << "\t" <<infantAllCauseMort() <<  lineEnd;
  }

  outputFile.close();

}

double Summary::infantAllCauseMort(){
  //DOCU

  double infantProductLimit=1.0;
  double* intervalSurviv = new double[Global::intervalsPerYear];
  for (size_t i=0;i<Global::intervalsPerYear; i++) {
    intervalSurviv[i]=(Global::infantIntervalsAtRisk[i]-Global::infantDeaths[i]);
    intervalSurviv[i]=intervalSurviv[i]/(Global::infantIntervalsAtRisk[i]);
    infantProductLimit=infantProductLimit*intervalSurviv[i];
  }
  delete[] intervalSurviv;
  return (1-infantProductLimit)*1000+_nonMalariaMortality;
}
void Summary::addToHost(double age, int value){

  int agegrp = this->ageGroup(age);
  _numHosts[_surveyPeriod-1][agegrp - 1]+=value;

}

void Summary::addToInfectedHost(double age, int value){

  int agegrp = this->ageGroup(age);
  _numInfectedHosts[_surveyPeriod-1][agegrp - 1]+=value;

}

void Summary::addToTotalInfections(double age, int value){

  int agegrp = this->ageGroup(age);
  _totalInfections[_surveyPeriod-1][agegrp - 1]+=value;
  
}


void Summary::addToTotalPatentInfections(double age, int value){

  int agegrp = this->ageGroup(age);
  _totalPatentInfections[_surveyPeriod-1][agegrp - 1]+=value;
}        

void Summary::addToPatentHost(double age, int value){

  int agegrp = this->ageGroup(age);
  _numPatentHosts[_surveyPeriod-1][agegrp - 1]+=value;
 
}
        
void Summary::addToSumLogDensity(double age, double value){

  int agegrp = this->ageGroup(age);
  _sumLogDensity[_surveyPeriod-1][agegrp - 1]+=value;

}    


void Summary::addToExpectedInfected(double age, double value){

  int agegrp = this->ageGroup(age);
  _numExpectedInfected[_surveyPeriod-1][agegrp - 1]+=value;

}

void Summary::addToPyrogenicThreshold(double age, double value){

  int agegrp = this->ageGroup(age);
  _pyrogenicThreshold[_surveyPeriod-1][agegrp - 1]+=value;

}
    
void Summary::addToSumX(double age, double value){

  int agegrp = this->ageGroup(age);
  _sumX[_surveyPeriod-1][agegrp - 1]+=value;

} 
  
void Summary::setAnnualAverageKappa(double kappa){

  _annualAverageKappa[_surveyPeriod - 1]=kappa; 

}

void Summary::setNumTransmittingHosts(double value){

  _numTransmittingHosts[_surveyPeriod -1]=value;

}

template <class T>
void writeArray(ostream& file, int measure, 
                int assimilatorMode, vector<T>& array){

  for (int i=0; i<(int)array.size()-1; i++){
 
      if (assimilatorMode ==  1)
        file << "\t" << array[i] << lineEnd;
      else
        file << i+1 << "\t" << 0 << "\t" 
             << measure << "\t" <<  array[i] <<  lineEnd;
      
    }

}

template <class T>
void writeArray(ostream& file, int measure, 
                int assimilatorMode, vector< vector<T> >& array){
  /*
    Write one of the summary arrays to the output file.
    unit: I/O unit
    measure: id of the measure recorded in this array (eg. nHosts,nPatent)
    sArray: summary array to be written to disk
  */ 

  for (int i=0; i< (int) array.size()-1; i++){
    for(int j=0; j< (int) array[i].size()-1; j++){
 
      if (assimilatorMode ==  1)
        file << "\t" << array[i][j] << lineEnd;
      else
        file << i+1 << "\t" << j+1 << "\t" 
             << measure << "\t" <<  array[i][j] <<  lineEnd;

    }
  }
  
}



