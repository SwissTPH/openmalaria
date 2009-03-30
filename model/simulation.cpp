/*
  This file is part of OpenMalaria.
 
  Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "simulation.h"

#include "boincWrapper.h"

#include "inputData.h"
#include "GSLWrapper.h"
#include "transmissionModel.h"
#include "population.h"
#include "summary.h"
#include "intervention.h"
#include <stdio.h>
#include "global.h"
#include "DescriptiveInfection.h"
/*

  Defines the duration of the simulation runs and the pre-timestep call-sequence for both the
  warmup and the main simulation.

  Duration of the main simulation (post-warmup), in 5 day intervals.
*/
int simulationDuration;
double relTimeInMainSim;

int Simulation::simulationTime;
Summary* Simulation::gMainSummary;
int Simulation::timeStep;

Simulation::Simulation(){
  //initialize input variables and allocate memory, correct order is important.
  validateInput();		// only checks get_mode() currently
  gMainSummary = new Summary();	// currently won't affect anything else
  Global::initGlobal();
  //FIXME: when multiple types exist, should happen through Infections
  DescriptiveInfection::initParameters();

  _transmissionModel = TransmissionModel::createTransmissionModel();
  
  gMainSummary->initSummaryParameters();
  _population = new Population(_transmissionModel, get_populationsize());
  EntoInterventionITN::initParameters();
  EntoInterventionIRS::initParameters();
  Human::initHumanParameters();
  IPTIntervention::initParameters();
  Vaccine::initParameters();
  Global::simulationMode=equilibriumMode;
  simulationDuration=get_simulation_duration();
  relTimeInMainSim=simulationDuration/(1.0*simulationDuration+Global::maxAgeIntervals);
  if (Global::modelVersion & INCLUDES_PK_PD) {
    initProteomeModule();
    initDrugModule(5*24*60, 24*60);
  }

}

Simulation::~Simulation(){
  //free memory
  gMainSummary->clearSummaryParameters();
  IPTIntervention::clearParameters();
  Vaccine::clearParameters();
  delete gMainSummary;
  delete _population;
  delete _transmissionModel;  
}

int Simulation::start(){
  int whichSeedfile;
  if (isCheckpoint()) {
    readCheckpoint();
  }
  else {
    simulationTime=0;
    _population->estimateRemovalRates();
    _population->initialiseHumanList();
  }
  _population->setupPyramid(isCheckpoint());
  if (isCheckpoint()){
    fstream checkpointFile;
    checkpointFile.open(checkpoint_name);
    checkpointFile >> whichSeedfile;
    checkpointFile.close();
    load_rng_state(whichSeedfile);
  }
  updateOneLifespan();
  _transmissionModel->initMainSimulation(get_populationsize());

  _population->initialiseInfantArrays();
  mainSimulation();
  return 0;
}

void Simulation::mainSimulation(){
  //TODO5D
  timeStep=0;
  gMainSummary->initialiseSummaries();
  Global::simulationMode=get_mode();
  
  while( timeStep <=  simulationDuration) {
    _population->implementIntervention(timeStep);
    //Calculate the current progress
    BoincWrapper::reportProgress(relTimeInMainSim*(timeStep/simulationDuration)+(1-relTimeInMainSim));
    //Here would be another place to write checkpoints. But then we need to save state of the surveys/events.
    ++simulationTime;
    ++timeStep;
    _population->update1();
    if ( timeStep ==  gMainSummary->getSurveyTimeInterval(gMainSummary->getSurveyPeriod()-1)) {
      _population->newSurvey();
    }
  }
  _population->clear();
  gMainSummary->writeSummaryArrays();
}

void Simulation::updateOneLifespan () {
  while( simulationTime <  Global::maxAgeIntervals) {
    BoincWrapper::reportProgress (simulationTime / Global::maxAgeIntervals * (1-relTimeInMainSim));
    if (BoincWrapper::timeToCheckpoint()) {
      writeCheckpoint();
      BoincWrapper::checkpointCompleted();
    }
    
    ++simulationTime;
    _population->update1();
  }
}


void Simulation::validateInput(){
  //all input validation goes here
}


bool Simulation::isCheckpoint(){
  //open the file in which the checkpoint would have been written
  fstream checkPonitFile(checkpoint_name,ios::in);

  if (!checkPonitFile.is_open()){
    //If File cannot be accessed:
    return false;
  }
  return true;

}

void Simulation::writeCheckpoint(){
  //set this to 1, so that the first checkpoint is written to file "checkpoint0"
  int lastCompleteCheckpoint = 1;

  //check if there already is a checkpoint. 
  fstream checkpointFile;
  checkpointFile.open(checkpoint_name, fstream::in);
  if(checkpointFile.is_open()){
    checkpointFile >> lastCompleteCheckpoint;
    checkpointFile.close();
  }
 
  /*
    Open the file in which the checkpoint data will be written
    use the one that is outdated
  */
  if ( lastCompleteCheckpoint ==  1) {
    f_checkpoint.open(checkpoint_0_name, ios::out | ios::binary);
    save_rng_state(0);
  }
  else {
    f_checkpoint.open(checkpoint_1_name, ios::out | ios::binary);
    save_rng_state(1);
  }
   
  f_checkpoint.precision(20);
  f_checkpoint << simulationTime << endl ;
  if (Global::modelVersion & INCLUDES_PK_PD) {
    f_checkpoint << ProteomeManager::getManager();
  }
  _population->writeLists(f_checkpoint);
  f_checkpoint.close();

  /*
    Now update the checkpoint file, to indicate that now the other datafile is
    up-to-date
  */
  fstream checkpointFile2;
  checkpointFile2.open(checkpoint_name,ios::out);
  if ( lastCompleteCheckpoint ==  1) {
    checkpointFile2 << "0" <<endl;
  }
  else {
    checkpointFile2 << "1" <<endl;
  }
  checkpointFile2.close();
}

void Simulation::readCheckpoint(){
    int lastCompleteCheckpoint;
    //Open the checkpoint file to see which cp datafile is up-to-date
    fstream checkpointFile;
    checkpointFile.open(checkpoint_name, ios::in);
    checkpointFile >> lastCompleteCheckpoint;
    checkpointFile.close();
    //Open the up-to-date file
    if ( lastCompleteCheckpoint ==  1){
      f_checkpoint.open(checkpoint_1_name, ios::in | ios::binary);
    }
    else {
      f_checkpoint.open(checkpoint_0_name, ios::in | ios::binary);
    }
    f_checkpoint >> simulationTime;
    if (Global::modelVersion & INCLUDES_PK_PD) {
      ProteomeManager* manager = ProteomeManager::getManager();
      f_checkpoint >> *manager;
    }
    _population->readLists(f_checkpoint);
    f_checkpoint.close();
}
