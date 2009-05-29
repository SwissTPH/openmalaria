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
/** Per-human interventions for stopping mosquitoes. */
#ifndef Hmod_EntoIntervention
#define Hmod_EntoIntervention

#include "global.h"

/** Base for anti-mosq. interventions. */
class EntoIntervention
{
public:
  /** Types of insecticide. */
  enum {
    NONE = 0,
  };
  
  EntoIntervention () : insecticide(NONE) {}
  
  friend ostream& operator<<(ostream& out, const EntoIntervention& ei);
  friend istream& operator>>(istream& in, EntoIntervention& ei);
  
  virtual void write(ostream& out) const;
  virtual void read(istream& in);
  
  /** Multiplier for host's availability to mosquitos. */
  virtual double availability() const =0;
  
protected:
  /** Gives a multiplier in the range [0,1] describing how
   * effective the intervention still is depending on its age.
   * 
   * Age is determined to be Simulation::simulationTime - dateOfUse. */
  virtual double decay() const =0;
  
  /** simulationTime - dateOfUse is the age of the intervention.
   * 
   * This is the date of last use. */
  int dateOfUse;
  
  /** Insecticide used. */
  int insecticide;
};

/** Insecticide Treated Nets (or untreated nets). */
class EntoInterventionITN : public EntoIntervention
{
public:
  /** Set some static values from XML. */
  static void initParameters();
  
  EntoInterventionITN () : netEffectiveness(0.0) {}
  
  void write(ostream& out) const;
  void read(istream& in);
  
  double availability() const;
  
  /** Multiplies the chance of a mosquito biting the host. */
  double probMosqBiting() const;
  
  /** Multiplies the probability of a mosquito finding a resting site after
   * biting. */
  double probMosqFindRestSite() const;
  
protected:
  double decay() const;
  
  // Alternately, could store brand/tdateOfUseype.
  
  /** Effectiveness of material */
  //NOTE: maybe unnecessary and should always be 1?
  double netEffectiveness;
  
  /** Life-span of net.
   * 
   * Days until it's totally useless.
   * 
   * Effectiveness tails off with age how? linearly? exponentially? */
  int netLifespan;
  
  /** The Weibull CDF is used to model decay. These are its constants: 1/λ and
   * k. */
  static double weibullDecayLambdaInv,
                weibullDecayk;
};

/** Indoor Residual Spraying */
class EntoInterventionIRS : public EntoIntervention
{
public:
  /** Set some static values from XML. */
  static void initParameters();
  
  double availability() const;
  
  /** Multiplies the chance of a mosquito resting succesfully. */
  double probMosqSurvivalResting() const;
  
private:
  double decay() const;
  
  /** Decay constant (see decay()). */
  static double decayLambdaInv;
};

ostream& operator<<(ostream& out, const EntoIntervention& ei);
istream& operator>>(istream& in, EntoIntervention& ei);

#endif
