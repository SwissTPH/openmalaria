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

#include "EntoIntervention.h"
#include <cmath>
#include "simulation.h"

double EntoInterventionITN::weibullDecayLambdaInv,
       EntoInterventionITN::weibullDecayk;
double EntoInterventionIRS::decayLambdaInv;

void EntoInterventionITN::initParameters () {
  //TODO: set from XML
  // Rough λ, k values to fit Polyester / Polyethylene: 7, 2.2 / 3.2, 1.8
  weibullDecayLambdaInv = 1.0/3.2;
  weibullDecayk = 1.8;
}
void EntoInterventionIRS::initParameters () {
  //TODO: set from XML
  decayLambdaInv = 1.;
}


/* ----- Main functionality ----- */

double EntoInterventionITN::decay () const {
  int age = Simulation::simulationTime - dateOfUse;
  
  return exp(-pow((double)age * weibullDecayLambdaInv, weibullDecayk));
}

double EntoInterventionITN::availability() const {
  // FIXME: number depends on net
  double E = 0.0;
  return 1.0 - (1.0 - E)*decay();
}

double EntoInterventionITN::probMosqBiting() const {
  // FIXME: number depends on net
  double E = 0.0;
  return 1.0 - (1.0 - E)*decay();
}

double EntoInterventionIRS::decay () const {
  int age = Simulation::simulationTime - dateOfUse;
  
  return exp (-age * decayLambdaInv);
}

double EntoInterventionIRS::availability() const {
  //FIXME: number depends on insecticide
  double scareConst = 0;	// probability of scaring mosquito off
  return 1.0 - scareConst * decay();
}

double EntoInterventionIRS::probMosqSurvivalResting() const {
  //FIXME: number depends on insecticide
  double killConst = 0;		// probability of killing mosquito
  return 1.0 - killConst * decay();
}


/* ----- Checkpointing code ----- */

void EntoInterventionITN::write(ostream& out) const {
  out << dateOfUse << endl;
  out << insecticide << endl;
  out << netEffectiveness << endl;
  out << netLifespan << endl;
}
void EntoInterventionITN::read(istream& in) {
  in >> dateOfUse;
  in >> insecticide;
  in >> netEffectiveness;
  in >> netLifespan;
}

void EntoIntervention::write(ostream& out) const {
  out << dateOfUse << endl;
  out << insecticide << endl;
}
void EntoIntervention::read(istream& in) {
  in >> dateOfUse;
  in >> insecticide;
}

ostream& operator<<(ostream& out, const EntoIntervention& ei) {
  //Not used yet
  //ei.write (out);
  return out;
}

istream& operator>>(istream& in, EntoIntervention& ei) {
  //ei.read (in);
  return in;
}
