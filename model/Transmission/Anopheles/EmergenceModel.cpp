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

#include "Transmission/Anopheles/EmergenceModel.h"
#include "Transmission/Anopheles/MosqTransmission.h"

#include <cmath>

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;

// -----  Initialisation of model, done before human warmup  ------

EmergenceModel::EmergenceModel() :
            EIRRotateAngle(numeric_limits<double>::quiet_NaN()),
            FSRotateAngle(numeric_limits<double>::quiet_NaN()),
            initNvFromSv(numeric_limits<double>::quiet_NaN()),
            initOvFromSv(numeric_limits<double>::quiet_NaN()),
            emergenceSurvival(1.0)
{
    forcedS_v.resize (SimTime::oneYear());
}

void EmergenceModel::initEIR(
    const scnXml::AnophelesParams& anoph,
    vector<double>& initialisationEIR,
    SimTime EIPDuration
){
    const scnXml::Seasonality& seasonality = anoph.getSeasonality();
    if ( seasonality.getInput() != "EIR" ) {
        throw util::xml_scenario_error("entomology.anopheles.seasonality.input: must be EIR (for now)");
        //TODO
    }
    // EIR for this species, with index 0 refering to value over first interval
    vecDay<double> speciesEIR (SimTime::oneYear());

    if ( seasonality.getFourierSeries().present() ) {
        const scnXml::FourierSeries& seasFC = seasonality.getFourierSeries().get();
        const scnXml::FourierSeries::CoefficSequence& fsCoeffic = seasFC.getCoeffic();

        FSCoeffic.reserve (2*fsCoeffic.size() + 1);

        FSCoeffic.push_back( 0.0 );     // value doesn't matter; EIR will be scaled
        for( scnXml::FourierSeries::CoefficConstIterator it=fsCoeffic.begin(); it!=fsCoeffic.end(); ++it ) {
            FSCoeffic.push_back( it->getA() );
            FSCoeffic.push_back( it->getB() );
        }
        // According to spec, EIR for first day of year (rather than EIR at the
        // exact start of the year) is generated with t=0 in Fourier series.
        EIRRotateAngle = seasFC.getEIRRotateAngle();
    } else if ( seasonality.getMonthlyValues().present() ) {
        const scnXml::MonthlyValues& seasM = seasonality.getMonthlyValues().get();
        if ( seasM.getSmoothing() != "fourier" ) {
            throw util::xml_scenario_error("entomology.anopheles.seasonality.monthlyValues.smoothing: only fourier supported at the moment");
            //TODO: should be easy to add no smoothing
        }

        const size_t N_m = 12;
        const scnXml::MonthlyValues::ValueSequence seq = seasM.getValue();
        assert( seq.size() == N_m );    // enforced by schema
        vector<double> months(N_m);
        double sum = 0.0;
        for( size_t i = 0; i < N_m; ++i ) {
            months[i] = seq[i];
            sum += months[i];
        }
        // arbitrary minimum we allow (cannot have zeros since we take the logarithm)
        double min = sum/1000.0;
        for( size_t i = 0; i < N_m; ++i ) {
            if ( months[i] < min )
                months[i] = min;
        }

        FSCoeffic.assign( 5, 0.0 );
        //TODO: determine whether to use Fourier Series Coefficient method or
        // Discrete Fourier Transform. Former is designed for integrals (which
        // roughly what we have?), the latter for discrete values. DFT doesn't
        // work well when number of intervals changes?
        vectors::logFourierCoefficients(months, FSCoeffic);

        // The above places the value for the first month at angle 0, so
        // effectively the first month starts at angle -2*pi/24 radians.
        // The value for the first day of the year should start 2*pi/(365*2)
        // radians later, so adjust EIRRotateAngle to compensate.
        EIRRotateAngle = M_PI * ( 1.0/12.0 - 1.0/365.0 );
    } else {
        assert( seasonality.getDailyValues().present() );      // XML loading code should enforce this
        throw util::xml_scenario_error("entomology.anopheles.seasonality.dailyValues: not supported yet");
        //TODO: can we encode as a Fourier series without smoothing?
    }

    if ( !seasonality.getAnnualEIR().present() ) {
        //TODO: work out when this is not required and implement code
        throw util::xml_scenario_error("entomology.anopheles.seasonality.annualEIR is required at the moment");
    }
    double targetEIR = seasonality.getAnnualEIR().get();
    
    // Now we rescale to get an EIR of targetEIR.
    // Calculate current sum as is usually done.
    vectors::expIDFT (speciesEIR, FSCoeffic, EIRRotateAngle);
    // And scale (also acts as a unit conversion):
    FSCoeffic[0] += log( targetEIR / vectors::sum( speciesEIR ) );
    
    // Calculate forced EIR for pre-intervention phase from FSCoeffic:
    vectors::expIDFT (speciesEIR, FSCoeffic, EIRRotateAngle);
    
    // Add to the TransmissionModel's EIR, used for the initalization phase.
    // Note: sum stays the same, units changes to per-time-step.
    for( SimTime i = SimTime::zero(); i < SimTime::oneYear(); i += SimTime::oneDay() ){
        // index 0 of initialisationEIR corresponds to first period of year
        initialisationEIR[mod_nn(i.inSteps(), SimTime::stepsPerYear())] += speciesEIR[i];
    }
    
#ifdef WITHOUT_BOINC
    if ( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ) {
        cout << "Annual EIR for "<<anoph.getMosquito()
             << ": "<<vectors::sum( speciesEIR )<<endl;
    }
#endif

    // Set other data used for mosqEmergeRate calculation:
    FSRotateAngle = EIRRotateAngle - (EIPDuration.inDays()+10)/365.*2.*M_PI;       // usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
    initNvFromSv = 1.0 / anoph.getPropInfectious();
    initOvFromSv = initNvFromSv * anoph.getPropInfected();
}

void EmergenceModel::scaleEIR( double factor ) {
    FSCoeffic[0] += log( factor );
}


// Every SimTime::oneTS() days:
void EmergenceModel::update () {
    emergenceSurvival = 1.0;
    for( size_t i = 0; i < emergenceReduction.size(); ++i ){
        emergenceSurvival *= 1.0 - emergenceReduction[i].current_value( sim::ts0() );
    }
}

void EmergenceModel::checkpoint (istream& stream){ (*this) & stream; }
void EmergenceModel::checkpoint (ostream& stream){ (*this) & stream; }

// -----  Summary and intervention functions  -----

void EmergenceModel::initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance ){
    if( emergenceReduction.size() <= instance )
        emergenceReduction.resize( instance+1 );
    
    if( elt.getEmergenceReduction().present() ){
        const scnXml::EmergenceReduction& elt2 = elt.getEmergenceReduction().get();
        if( elt2.getInitial() > 1.0 )
            throw util::xml_scenario_error( "emergenceReduction intervention: initial effect must be ≤ 1" );
        emergenceReduction[instance].set (elt2.getInitial(), elt2.getDecay(), "emergenceReduction");
    }
}

void EmergenceModel::deployVectorPopInterv (size_t instance) {
    assert( instance < emergenceReduction.size() );
    emergenceReduction[instance].deploy( sim::now() );
}


}
}
}
