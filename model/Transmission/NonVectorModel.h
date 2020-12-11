/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * 
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef Hmod_NonVectorModel
#define Hmod_NonVectorModel

#include "TransmissionModel.h"

namespace scnXml {
  class NonVector;
}

namespace OM { namespace Transmission {
    
//! Base transmission model, as used in Phase A
class NonVectorModel : public TransmissionModel { 
public:

  NonVectorModel(const scnXml::Entomology& entoData, const scnXml::NonVector& nonVectorData);
  virtual ~NonVectorModel();
  
  virtual void init2 (const Population& population);
  
  virtual void initVectorInterv( const scnXml::Description::AnophelesSequence& list,
        size_t instance, const string& name );
  virtual void initVectorTrap( const scnXml::VectorTrap::DescriptionSequence list,
        size_t instance, const scnXml::VectorTrap::NameOptional name );
  virtual void initNonHumanHostsInterv( const scnXml::Description2::AnophelesSequence list,
        const scnXml::DecayFunction& decay, size_t instance, const string& name );
  virtual void initAddNonHumanHostsInterv( const scnXml::Description3::AnophelesSequence list,
        const string& name );

  virtual void scaleEIR (double factor);
//   virtual void scaleXML_EIR (scnXml::Entomology&, double factor) const;
  
  virtual SimTime minPreinitDuration ();
  virtual SimTime expectedInitDuration ();
  virtual SimTime initIterate ();
  
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
  virtual void changeEIRIntervention (const scnXml::NonVector&);
  
  virtual void deployVectorPopInterv (size_t instance);
  virtual void deployVectorTrap( size_t instance, double popSize, SimTime lifespan );
  virtual void deployNonHumanHostsInterv( size_t instance, string name );
  virtual void deployAddNonHumanHosts( string name, double popSize, SimTime lifespan );

  virtual void uninfectVectors();
  
  virtual void vectorUpdate (const Population& population) {}
  virtual void update (const Population& population);
  virtual void calculateEIR(OM::Host::Human& human, double ageYears, vector< double >& EIR) const;
  
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
  SimTime nSpore;
  //@}
  
  /** EIR per time interval during the intervention period. Value at index
   * sim::intervTime().inSteps() used each time-step.
   * 
   * Units: inoculations per adult per time step */
  vector<double> interventionEIR;
  
  /** When simulationMode == dynamicEIR, this is the annual cycle of kappa
   * from the warmup phase and has length 1 year (in time steps).
   * 
   * When simulationMode == equilibriumMode, this may be multiple years long and
   * is used to collect values of kappa (human infectiousness).
   * 
   * In either case, sim::ts0().moduloSteps(initialKappa.size()) is the index
   * for the current infectiousness during updates. */
  vector<double> initialKappa; 
};
} }
#endif
