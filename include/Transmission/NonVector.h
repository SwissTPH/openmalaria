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
#ifndef Hmod_NonVectorTransmission
#define Hmod_NonVectorTransmission

#include "TransmissionModel.h"

namespace scnXml {
  class NonVector;
}

namespace OM { namespace Transmission {
    
//! Base transmission model, as used in Phase A
class NonVectorTransmission : public TransmissionModel { 
public:

  NonVectorTransmission(const scnXml::NonVector& nonVectorData);
  virtual ~NonVectorTransmission();
  
  virtual void scaleEIR (double factor);
  virtual void scaleXML_EIR (scnXml::EntoData&, double factor) const;
  
  virtual TimeStep minPreinitDuration ();
  virtual TimeStep expectedInitDuration ();
  virtual TimeStep initIterate ();
  
  /** Change the scnXml::NonVector data (changeEIR intervention). */
  void setNonVectorData (const scnXml::NonVector& nonVectorData);
  
  /** @brief Sets up the EIR used in a change of EIR intervention.
   * 
   * EIR is always set from intervention-period step 0, not the current step.
   *
   * Reads in the estimates of the EIR for each village and each day
   * and converts this into EIR estimates per five day period
   * assuming that the annual cycle repeated during the pre-intervention period
   * 
   * Similar calculation to that used during initialization. */
  void setTransientEIR (const scnXml::NonVector& nonVectorData); 
  
  virtual void changeEIRIntervention (const scnXml::NonVector&);
  
  virtual void uninfectVectors();
  
  virtual void update (const std::list<Host::Human>& population, int populationSize);
  virtual double calculateEIR(PerHostTransmission& perHost, double ageYears);
  
private:

  /// Processes each daily EIR estimate, allocating each day in turn to the
  /// appropriate time period. EIRdaily is the value of the daily EIR read in
  /// from the .XML file.
  void updateEIR (int day, double EIRdaily); 
  double averageEIR (const scnXml::NonVector& nonVectorData); 
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
//! multiplier used to calculate a positive EIR value where the measured value is zero
/* 
  0.01 was old pv(30) Now a constant. min_EIR_mult multiplies the average EIR to obtain a value used for the EIR during periods when it is too low 
  to be measureable. The value of 0.01 was old pv(30) Now a constant. 
  0.01 was old pv(30) Now a constant. 
*/ 
  static const double min_EIR_mult; 

  /** @brief Variables set by constructor.
   *
   * There shouldn't be any need to checkpoint these, at least before
   * interventions take effect. */
  //@{
  //! Variance of Infection Rate according to fielddata 
  static const double totalInfectionrateVariance; 
  
  //! The duration of sporogony in time steps
  // doesn't need checkpointing
  TimeStep nspore;
  //@}
  
  /// EIR per time interval during the intervention period
  /// Units: inoculations per adult per timestep
  vector<double> interventionEIR;
  
  /** When simulationMode == dynamicEIR, this is the annual cycle of kappa
   * from the warmup phase and has length 1 year (in timesteps). Index for this
   * time-step is TimeStep::simulation % initialKappa.size().
   * 
   * When simulationMode == equilibriumMode, this may be multiple years long and
   * is used to collect values of kappa (human infectiousness). */
  vector<double> initialKappa; 
};
} }
#endif
