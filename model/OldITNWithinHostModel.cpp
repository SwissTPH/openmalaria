/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#include "OldITNWithinHostModel.h"
#include "human.h"
#include "simulation.h"
#include "GSLWrapper.h"

void OldITNWithinHostModel::SPAction(Human& human){
  /*TODO if we want to look at presumptive SP treatment with the PkPD model we
  need to add some code here that will be conditionally implemented depending on the
  model version.*/

  double rnum;
  std::list<DescriptiveInfection>::iterator i=infections.begin();
  while(i != infections.end()){
    if ( 1+Simulation::simulationTime-i->getStartDate()-Global::latentp > 0){
      rnum=W_UNIFORM();
      if ((rnum<=IPTIntervention::genotypeACR[i->getGenoTypeID()-1]) &&
           (Simulation::simulationTime - human.getLastSPDose() <= IPTIntervention::genotypeProph[i->getGenoTypeID()-1])) {
        i->destroy();
        i=infections.erase(i);
        _MOI--;
           }
           else{
             i++;
           }
    }
    else{
      i++;
    }
  }
}

void OldITNWithinHostModel::IPTattenuateAsexualDensity (std::list<DescriptiveInfection>::iterator i) {
  if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY) {
    if (i->getSPattenuate() == 1) {
      i->multiplyDensity(1.0/IPTIntervention::genotypeAtten[i->getGenoTypeID() - 1]);
      timeStepMaxDensity=(double)timeStepMaxDensity/IPTIntervention::genotypeAtten[i->getGenoTypeID() - 1];
      _SPattenuationt=(int)std::max(_SPattenuationt*1.0, (i->getStartDate()+(i->getDuration()/Global::interval) * IPTIntervention::genotypeAtten[i->getGenoTypeID() - 1]));
    }
  }
}

void OldITNWithinHostModel::IPTattenuateAsexualMinTotalDensity (Human& human) {
  if (Global::modelVersion & ATTENUATION_ASEXUAL_DENSITY) {
    if (_SPattenuationt > Simulation::simulationTime &&  human.getTotalDensity() <  10) {
      human.setTotalDensity(10);
      human.setCumulativeY(human.getCumulativeY()+10);
    }
  }
}


// -----  Data checkpointing  -----

void OldITNWithinHostModel::read(istream& in) {
  readOWHM (in);
  in >> _SPattenuationt;
}
void OldITNWithinHostModel::write(ostream& out) const {
  writeOWHM (out);
  out << _SPattenuationt << endl;
}
