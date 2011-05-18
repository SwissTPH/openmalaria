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
#include "Transmission/PerHostTransmission.h"
#include "Transmission/Vector/VectorTransmission.h"	//TODO: remove dependency
#include "inputData.h"

namespace OM { namespace Transmission {
        using namespace OM::util;

// -----  PerHostTransmission static  -----

AgeGroupInterpolation* PerHostTransmission::relAvailAge = AgeGroupInterpolation::dummyObject();
shared_ptr<DecayFunction> PerHostTransmission::ITNDecay;
shared_ptr<DecayFunction> PerHostTransmission::IRSDecay;
shared_ptr<DecayFunction> PerHostTransmission::VADecay;

void PerHostTransmission::init () {
    relAvailAge = AgeGroupInterpolation::makeObject( InputData().getModel().getHuman().getAvailabilityToMosquitoes(), "availabilityToMosquitoes" );
    
    /*FIXME
    scnXml::Descriptions intervDesc = InputData().getInterventions().getDescriptions();
    if ( InputData.isInterventionActive(Interventions::ITN) ) {
        if( !intervDesc.getITNDecay().present() ){
            throw util::xml_scenario_error ("ITN intervention without description of decay");
        }
        ITNDecay = DecayFunction::makeObject( intervDesc.getITNDecay().get(), "ITNDecay" );
    }
    if ( InputData.isInterventionActive(Interventions::IRS) ) {
        if( !intervDesc.getIRSDecay().present() ){
            throw util::xml_scenario_error ("IRS intervention without description of decay");
        }
        IRSDecay = DecayFunction::makeObject( intervDesc.getIRSDecay().get(), "IRSDecay" );
    }
    if ( InputData.isInterventionActive(Interventions::VEC_AVAIL) ) {
        if( !intervDesc.getVADecay().present() ){
            throw util::xml_scenario_error ("Vector availability intervention without description of decay");
        }
        VADecay = DecayFunction::makeObject( intervDesc.getVADecay().get(), "VADecay" );
    }
    */
}
void PerHostTransmission::cleanup (){
    AgeGroupInterpolation::freeObject( relAvailAge );
}


// -----  PerHostTransmission non-static -----

PerHostTransmission::PerHostTransmission () :
    outsideTransmission(false),
    timestepITN(TimeStep::never),
    timestepIRS(TimeStep::never),
    timestepVA(TimeStep::never)
{
    if( ITNDecay.get() != 0 )
        hetSampleITN = ITNDecay->hetSample();
    if( IRSDecay.get() != 0 )
        hetSampleIRS = IRSDecay->hetSample();
    if( VADecay.get() != 0 )
        hetSampleVA = VADecay->hetSample();
}
void PerHostTransmission::initialise (TransmissionModel& tm, double availabilityFactor) {
  _relativeAvailabilityHet = availabilityFactor;
  VectorTransmission* vTM = dynamic_cast<VectorTransmission*> (&tm);
  if (vTM) {
    species.resize (vTM->numSpecies);
    for (size_t i = 0; i < vTM->numSpecies; ++i)
      species[i].initialise (vTM->species[i].getHumanBase(), availabilityFactor);
  }
}


// Note: in the case an intervention is not present, we can use the approximation
// of Weibull decay over (TimeStep::simulation - TimeStep::never) timesteps
// (easily large enough for conceivable Weibull params that the value is 0.0 when
// rounded to a double. Performance-wise it's perhaps slightly slower than using
// an if() when interventions aren't present.
double PerHostTransmission::entoAvailabilityHetVecItv (const HostCategoryAnopheles& base, size_t speciesIndex) const {
  double alpha_i = species[speciesIndex].entoAvailability;
  if (timestepITN >= TimeStep(0))
    alpha_i *= (1.0 - base.ITNDeterrency * ITNDecay->eval (TimeStep::simulation - timestepITN, hetSampleITN));
  if (timestepIRS >= TimeStep(0))
    alpha_i *= (1.0 - base.IRSDeterrency * IRSDecay->eval (TimeStep::simulation - timestepIRS, hetSampleIRS));
  if (timestepVA >= TimeStep(0))
    alpha_i *= (1.0 - base.VADeterrency * VADecay->eval (TimeStep::simulation - timestepVA, hetSampleVA));

  return alpha_i;
}
double PerHostTransmission::probMosqBiting (const HostCategoryAnopheles& base, size_t speciesIndex) const {
  double P_B_i = species[speciesIndex].probMosqBiting;
  if (timestepITN >= TimeStep(0))
    P_B_i *= (1.0 - base.ITNPreprandialKillingEffect * ITNDecay->eval (TimeStep::simulation - timestepITN, hetSampleITN));
  return P_B_i;
}
double PerHostTransmission::probMosqResting (const HostCategoryAnopheles& base, size_t speciesIndex) const {
  double P_C_i = species[speciesIndex].probMosqFindRestSite;
  if (timestepITN >= TimeStep(0))
    P_C_i *= (1.0 - base.ITNPostprandialKillingEffect * ITNDecay->eval (TimeStep::simulation - timestepITN, hetSampleITN));
  double P_D_i = species[speciesIndex].probMosqSurvivalResting;
  if (timestepIRS >= TimeStep(0))
    P_D_i *= (1.0 - base.IRSKillingEffect * IRSDecay->eval (TimeStep::simulation - timestepIRS, hetSampleIRS));
  return P_C_i * P_D_i;
}


// ----- HostMosquitoInteraction non-static -----

void HostMosquitoInteraction::initialise (HostCategoryAnopheles& base, double availabilityFactor)
{
  //TODO: vary to simulate heterogeneity

  entoAvailability = base.entoAvailability * availabilityFactor;
  probMosqBiting = base.probMosqBiting;
  probMosqFindRestSite = base.probMosqFindRestSite;
  probMosqSurvivalResting = base.probMosqSurvivalResting;
}

} }
