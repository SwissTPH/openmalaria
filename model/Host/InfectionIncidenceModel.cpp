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

#include "Host/InfectionIncidenceModel.h"
#include "Host/Human.h"
#include "Transmission/PerHost.h"
#include "WithinHost/WHInterface.h"
#include "Parameters.h"
#include "Monitoring/Continuous.h"
#include "util/ModelOptions.h"
#include "util/random.h"
#include "util/errors.h"

#include <stdexcept>
#include <cmath>

namespace OM { namespace Host {
    using namespace OM::util;

// ———  model constants set by init()  ———

/* Shape constant of (Gamma) distribution of availability
real, parameter :: BaselineAvailabilityGammaShapeParam =1.0 */
double baseline_avail_shape_param, baseline_avail_offset;

//VARIABLES INCLUDED IN CORE GETs of number of infections 
//! Describes the shape of the Infectionrate distribution, related to the baseline availabilty distr.
double inf_rate_shape_param, inf_rate_offset;

/** @brief Parameters used in the "expected number of infections" model */
//@{
//!Steepness of relationship between success of inoculation and Xp in Phase A model 
double gamma_p;
//!Lower limit of success probability of inoculations at high exposure in Phase A model 
double Sinf;
//!Lower limit of success probability of inoculations in immune individuals in Phase A model 
double Simm;
//!1 over the critical value of cumulative number of entomologic inoculations in Phase A model 
double Xstar_pInv;
//!1 over the critical value of EIR in Phase A pre-erythrocytic model 
double EstarInv;
//@}

// model options:
bool opt_neg_bin_mass_action = false, opt_lognormal_mass_action = false,
        opt_no_pre_erythrocytic = false, opt_any_het = false;

// ———  variables  ———
int InfectionIncidenceModel::ctsNewInfections = 0;

// -----  static initialisation  -----

void InfectionIncidenceModel::init ( const Parameters& parameters ) {
    baseline_avail_shape_param = parameters[Parameters::BASELINE_AVAILABILITY_SHAPE];
    baseline_avail_offset = 0.5 * pow(baseline_avail_shape_param, 2);
    
    gamma_p=parameters[Parameters::GAMMA_P];
    Sinf=1-exp(-parameters[Parameters::NEG_LOG_ONE_MINUS_SINF]);
    Simm=parameters[Parameters::SIMM];
    EstarInv = 1.0/parameters[Parameters::E_STAR];
    Xstar_pInv = 1.0/parameters[Parameters::X_STAR_P];
    
    //! constant defining the constraint for the Gamma shape parameters
    /// Used for the case where availability is assumed gamma distributed
    double r_square_Gamma;
    /*
    //! Expected number of inoculations
    /// Product of measured EIR, susceptibility and length of time Global::interval
    double gsi = 1.0;
    
    r_square_Gamma=(totalInfectionrateVariance**2-gsi*BaselineAvailabilityMean)/(gsi*BaselineAvailabilityMean)**2
    r_square_Gamma must be greater than zero, so r_square_LogNormal is also. 
    */
    r_square_Gamma=0.649; //such that r_square_LogNormal =0.5
    
    opt_no_pre_erythrocytic = util::ModelOptions::option (util::NO_PRE_ERYTHROCYTIC);
    opt_neg_bin_mass_action = util::ModelOptions::option (util::NEGATIVE_BINOMIAL_MASS_ACTION);
    if (opt_neg_bin_mass_action) {
        inf_rate_shape_param = (baseline_avail_shape_param+1.0) / (r_square_Gamma*baseline_avail_shape_param - 1.0);
        inf_rate_shape_param=std::max(inf_rate_shape_param, 0.0);
    } else if (util::ModelOptions::option (util::LOGNORMAL_MASS_ACTION)) {
        opt_lognormal_mass_action = true;
        
        //! constant defining the constraint for the log Normal variance
        /// Used for the case where availability is assumed log Normally distributed
        double r_square_LogNormal = log(1.0+r_square_Gamma);
        
        inf_rate_shape_param = sqrt(r_square_LogNormal - 1.86 * pow(baseline_avail_shape_param, 2));
        inf_rate_shape_param=std::max(inf_rate_shape_param, 0.0);
        inf_rate_offset = 0.5 * pow(inf_rate_shape_param, 2);
        if( (boost::math::isnan)(inf_rate_shape_param) ){
            throw util::xml_scenario_error( "bad parameter 16 (BASELINE_AVAILABILITY_SHAPE)" );
        }
    }else{
        if( util::ModelOptions::option( util::TRANS_HET ) ||
            util::ModelOptions::option( util::COMORB_TRANS_HET ) ||
            util::ModelOptions::option( util::TRANS_TREAT_HET ) ||
            util::ModelOptions::option( util::TRIPLE_HET ) )
        {
            opt_any_het = true;
            cerr << "Warning: will use heterogeneity workaround." << endl;
        }
    }
    
    Monitoring::Continuous.registerCallback( "new infections", "\tnew infections", &InfectionIncidenceModel::ctsReportNewInfections );
}


// -----  other static methods  -----

void InfectionIncidenceModel::ctsReportNewInfections (ostream& stream){
    stream << '\t' << ctsNewInfections;
    ctsNewInfections = 0;
}


// -----  non-static non-checkpointing constructors  -----

InfectionIncidenceModel* InfectionIncidenceModel::createModel () {
    if (opt_neg_bin_mass_action) return new NegBinomMAII ();
    if(opt_lognormal_mass_action) return new LogNormalMAII ();
    if (opt_any_het) return new HeterogeneityWorkaroundII();
    return new InfectionIncidenceModel();
}

InfectionIncidenceModel::InfectionIncidenceModel () :
  m_pInfected(0.0), m_cumulativeEIRa(0.0)
{}


// -----  non-static methods  -----

double InfectionIncidenceModel::getAvailabilityFactor(double baseAvailability) {
  return baseAvailability;
}
double NegBinomMAII::getAvailabilityFactor(double baseAvailability) {
    // Gamma sample with k=BaselineAvailabilityShapeParam; mean is baseAvailability
  return random::gamma(baseline_avail_shape_param,
		 baseAvailability/baseline_avail_shape_param);
}
double LogNormalMAII::getAvailabilityFactor(double baseAvailability) {
    // given BaselineAvailabilityShapeParam = sqrt (log (1 + variance/mean²))
    // and baseAvailability = mean, this is a draw from the log-normal distribution.
    if( baseAvailability != 1.0 ){
        // NOTE: shouldn't the normal_mean parameter be adjusted when baseAvailability != 1.0?
        throw TRACED_EXCEPTION_DEFAULT("LogNormalMAII::getAvailabilityFactor");
    }
  return random::log_normal (log(baseAvailability) - baseline_avail_offset,
		     baseline_avail_shape_param);
}

void InfectionIncidenceModel::summarize (const Host::Human& human) {
    mon::reportStatMHF( mon::MHF_EXPECTED_INFECTED, human, m_pInfected );
}


double InfectionIncidenceModel::getModelExpectedInfections (double effectiveEIR, const Transmission::PerHost& phTrans) {
  // First two lines are availability adjustment: S_1(i,t) from AJTMH 75 (suppl 2) p12 eqn. (5)
  // Note that NegBinomMAII and LogNormalMAII supercede this model; see below
  return (Sinf+(1-Sinf) / 
    (1 + effectiveEIR/SimTime::oneTS().inDays()*EstarInv)) *
    susceptibility() * effectiveEIR;
}
double HeterogeneityWorkaroundII::getModelExpectedInfections (double effectiveEIR, const Transmission::PerHost& phTrans) {
  return (Sinf+(1-Sinf) / 
    (1 + effectiveEIR/(SimTime::oneTS().inDays()*phTrans.relativeAvailabilityHet())*EstarInv)) *
    susceptibility() * effectiveEIR;
}
double NegBinomMAII::getModelExpectedInfections (double effectiveEIR, const Transmission::PerHost&) {
  // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pmed.1001157.s009
  return random::gamma(inf_rate_shape_param,
      effectiveEIR * susceptibility() / inf_rate_shape_param);
}
double LogNormalMAII::getModelExpectedInfections (double effectiveEIR, const Transmission::PerHost&) {
  // Documentation: http://www.plosmedicine.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pmed.1001157.s009
    //TODO: is this equivalent to gsl_ran_lognormal?
  return random::sampleFromLogNormal(random::uniform_01(),
      log(effectiveEIR * susceptibility()) - inf_rate_offset, inf_rate_shape_param);
}

double InfectionIncidenceModel::susceptibility () {
  if (opt_no_pre_erythrocytic) {
    //! The average proportion of bites from sporozoite positive mosquitoes resulting in infection. 
    /*! 
    This is computed as 0.19 (the value S from a neg bin mass action model fitted 
    to Saradidi data, divided by 0.302 (the ratio of body surface area in a 
    0.5-6 year old child (as per Saradidi) to adult) 
    \sa getExpectedNumberOfInfections() 
    */ 
    return 0.702;
  } else {
    // S_2(i,t) from AJTMH 75 (suppl 2) p12 eqn. (7)
    return Simm + (1.0-Simm) /
      (1.0 + pow(m_cumulativeEIRa*Xstar_pInv, gamma_p));
  }
}

int InfectionIncidenceModel::numNewInfections (const Human& human, double effectiveEIR) {
  double expectedNumInfections = getModelExpectedInfections (effectiveEIR, human.perHostTransmission);
  // error check (should be OK if kappa is checked, for nonVector model):
  if( !(boost::math::isfinite)(effectiveEIR) ){
    ostringstream out;
    out << "effectiveEIR is not finite: " << effectiveEIR << endl;
    throw TRACED_EXCEPTION (out.str(), util::Error::EffectiveEIR);
  }
  
  //Introduce the effect of vaccination. Note that this does not affect cumEIR.
    expectedNumInfections *= human.getVaccine().getFactor( interventions::Vaccine::PEV );
  
  //Update pre-erythrocytic immunity
  m_cumulativeEIRa+=effectiveEIR;
  
  m_pInfected = 1.0 - exp(-expectedNumInfections) * (1.0-m_pInfected);
  if (m_pInfected < 0.0)
    m_pInfected = 0.0;
  else if (m_pInfected > 1.0)
    m_pInfected = 1.0;
  
  if (expectedNumInfections > 0.0000001){
    int n = random::poisson(expectedNumInfections);
    if( n > WithinHost::WHInterface::MAX_INFECTIONS ){
        n = WithinHost::WHInterface::MAX_INFECTIONS;
    }
    mon::reportEventMHI( mon::MHR_NEW_INFECTIONS, human, n );
    ctsNewInfections += n;
    return n;
  }
  if ( (boost::math::isnan)(expectedNumInfections) ){	// check for not-a-number
      // bad Params::BASELINE_AVAILABILITY_SHAPE ?
      throw TRACED_EXCEPTION( "numNewInfections: NaN", util::Error::NumNewInfections );
  }
  return 0;
}

} }
