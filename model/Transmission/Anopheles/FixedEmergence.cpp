/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2012 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Anopheles/FixedEmergence.h"
#include "Transmission/Anopheles/Transmission.h"
#include "Transmission/Anopheles/Nv0DelayFitting.h"

#include "util/vectors.h"
#include "util/CommandLine.h"
#include "util/errors.h"

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;


// -----  Initialisation of model, done before human warmup  ------

void FixedEmergence::initEIR(
    const scnXml::AnophelesParams& anoph,
    vector<double>& initialisationEIR,
    int EIPDuration
){
    const scnXml::Seasonality& seasonality = anoph.getSeasonality();
    if ( seasonality.getInput() != "EIR" ) {
        throw util::xml_scenario_error("entomology.anopheles.seasonality.input: must be EIR (for now)");
        //TODO
    }
    // EIR for this species, with index 0 refering to value over first interval
    vector<double> speciesEIR (TimeStep::DAYS_IN_YEAR);

    if ( seasonality.getFourierSeries().present() ) {
        const scnXml::FourierSeries& seasFC = seasonality.getFourierSeries().get();
        const scnXml::FourierSeries::CoefficSequence& fsCoeffic = seasFC.getCoeffic();

        FSCoeffic.reserve (2*fsCoeffic.size() + 1);

        FSCoeffic.push_back( 0.0 );     // value doesn't matter; EIR will be scaled
        for ( scnXml::FourierSeries::CoefficConstIterator it=fsCoeffic.begin(); it!=fsCoeffic.end(); ++it ) {
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
        double months[N_m];
        double sum = 0.0;
        for ( size_t i = 0; i < N_m; ++i ) {
            months[i] = seq[i];
            sum += months[i];
        }
        // arbitrary minimum we allow (cannot have zeros since we take the logarithm)
        double min = sum/1000.0;
        for ( size_t i = 0; i < N_m; ++i ) {
            if ( months[i] < min )
                months[i] = min;
        }

        const double PI = 3.14159265;
        const double w = 2.0 * PI / N_m;
        FSCoeffic.assign( 5, 0.0 );

        // Note: we use our values as the left-hand-side of our regions
        for ( size_t i = 0; i < N_m; ++i ) {
            double val = log( months[i] );
            FSCoeffic[0] += val;
            FSCoeffic[1] += val * cos( w*i );
            FSCoeffic[2] += val * sin( w*i );
            FSCoeffic[3] += val * cos( 2.0*w*i );
            FSCoeffic[4] += val * sin( 2.0*w*i );
        }
        FSCoeffic[0] /=N_m;
        FSCoeffic[1] *= 2.0 / N_m;
        FSCoeffic[2] *= 2.0 / N_m;
        FSCoeffic[3] *= 2.0 / N_m;
        FSCoeffic[4] *= 2.0 / N_m;

        // The above places the value for the first month at angle 0, so
        // effectively the first month starts at angle -2*pi/24 radians.
        // The value for the first day of the year should start 2*pi/(365*2)
        // radians later, so adjust EIRRotateAngle to compensate.
        EIRRotateAngle = M_PI * ( 1.0/12.0 - 1.0/365.0 );
    } else {
        assert( seasonality.getDailyValues().present() );      // XML loading code should enforce this
        throw util::xml_scenario_error("entomology.anopheles.seasonality.dailyValues: not supported yet");
        //TODO
    }

    if ( !seasonality.getAnnualEIR().present() ) {
        //TODO: work out when this is not required and implement code
        throw util::xml_scenario_error("entomology.anopheles.seasonality.annualEIR is required at the moment");
    }
    double targetEIR = seasonality.getAnnualEIR().get();

    // Now we rescale to get an EIR of targetEIR.
    // Calculate current sum as is usually done.
    vectors::calcExpFourierSeries (speciesEIR, FSCoeffic, EIRRotateAngle);
    // And scale (also acts as a unit conversion):
    FSCoeffic[0] += log( targetEIR / vectors::sum( speciesEIR ) );

    // Calculate forced EIR for pre-intervention phase from FSCoeffic:
    vectors::calcExpFourierSeries (speciesEIR, FSCoeffic, EIRRotateAngle);

    // Add to the TransmissionModel's EIR, used for the initalization phase:
    for (int i = 0; i < TimeStep::DAYS_IN_YEAR; ++i) {
        // index 1 of initialisationEIR corresponds to first period of year
        initialisationEIR[(1 + i / TimeStep::interval) % TimeStep::stepsPerYear] += speciesEIR[i];
    }

    if ( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ) {
        cout << "Annual EIR for "<<anoph.getMosquito()
             << ": "<<vectors::sum( speciesEIR )<<endl;
    }

    // Set other data used for mosqEmergeRate calculation:
    FSRotateAngle = EIRRotateAngle - (EIPDuration+10)/365.*2.*M_PI;       // usually around 20 days; no real analysis for effect of changing EIPDuration or mosqRestDuration
    initNvFromSv = 1.0 / anoph.getPropInfectious();
    initNv0FromSv = initNvFromSv * anoph.getPropInfected();       // temporarily use of initNv0FromSv
    
    // -----  allocate memory  -----
    quinquennialS_v.assign (TimeStep::fromYears(5).inDays(), 0.0);
    forcedS_v.resize (TimeStep::DAYS_IN_YEAR);
    mosqEmergeRate.resize (TimeStep::DAYS_IN_YEAR); // Only needs to be done here if loading from checkpoint
}

void FixedEmergence::scaleEIR( double factor ) {
    FSCoeffic[0] += log( factor );
}


// -----  Initialisation of model which is done after creating initial humans  -----

void FixedEmergence::init2( double tsP_A, double tsP_df, double EIRtoS_v, Transmission& transmission ){
    // -----  Calculate required S_v based on desired EIR  -----
    
    double initOvFromSv = initNv0FromSv;  // temporarily use of initNv0FromSv
    initNv0FromSv = initNvFromSv * (1.0 - tsP_A - tsP_df);

    // We scale FSCoeffic to give us S_v instead of EIR.
    // Log-values: adding log is same as exponentiating, multiplying and taking
    // the log again.
    FSCoeffic[0] += log( EIRtoS_v);
    vectors::calcExpFourierSeries (forcedS_v, FSCoeffic, FSRotateAngle);
    
    transmission.initState ( tsP_A, tsP_df, initNvFromSv, initOvFromSv, forcedS_v );
    
    // Crude estimate of mosqEmergeRate: (1 - P_A(t) - P_df(t)) / (T * ρ_S) * S_T(t)
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);
    
    // All set up to drive simulation from forcedS_v
}


// -----  Initialisation of model which is done after running the human warmup  -----

bool FixedEmergence::initIterate (Transmission& transmission) {
    // Try to match S_v against its predicted value. Don't try with N_v or O_v
    // because the predictions will change - would be chasing a moving target!
    // EIR comes directly from S_v, so should fit after we're done.

    double factor = vectors::sum (forcedS_v)*5 / vectors::sum(quinquennialS_v);
    //cout << "Pre-calced Sv, dynamic Sv:\t"<<sumAnnualForcedS_v<<'\t'<<vectors::sum(annualS_v)<<endl;
    if (!(factor > 1e-6 && factor < 1e6)) {
        if ( vectors::sum(forcedS_v) == 0.0 ) {
            return false;   // no EIR desired: nothing to do
        }
        cerr << "Input S_v for this vector:\t"<<vectors::sum(forcedS_v)<<endl;
        cerr << "Simulated S_v:\t\t\t"<<vectors::sum(quinquennialS_v)/5.0<<endl;
        throw TRACED_EXCEPTION ("factor out of bounds (likely a code error)",util::Error::VectorFitting);
    }

    //cout << "Vector iteration: adjusting with factor "<<factor<<endl;
    // Adjusting mosqEmergeRate is the important bit. The rest should just
    // bring things to a stable state quicker.
    initNv0FromSv *= factor;
    initNvFromSv *= factor;     //(not currently used)
    vectors::scale (mosqEmergeRate, factor);
    transmission.initIterateScale (factor);
    vectors::scale (quinquennialS_v, factor); // scale so we can fit rotation offset

    // average annual period of S_v over 5 years
    vector<double> avgAnnualS_v( TimeStep::fromYears(1).inDays(), 0.0 );
    for ( int i = 0; i < TimeStep::fromYears(5).inDays(); ++i ) {
        avgAnnualS_v[i % TimeStep::fromYears(1).inDays()] =
            quinquennialS_v[i] / 5.0;
    }

    // Once the amplitude is approximately correct, we try to find a
    // rotation offset.
    double rAngle = Nv0DelayFitting::fit<double> (EIRRotateAngle, FSCoeffic, avgAnnualS_v);
    //cout << "Vector iteration: rotating with angle (in radians): " << rAngle << endl;
    // annualS_v was already rotated by old value of FSRotateAngle, so increment:
    FSRotateAngle -= rAngle;
    vectors::calcExpFourierSeries (forcedS_v, FSCoeffic, FSRotateAngle);
    // We use the stored initXxFromYy calculated from the ideal population age-structure (at init).
    mosqEmergeRate = forcedS_v;
    vectors::scale (mosqEmergeRate, initNv0FromSv);

    const double LIMIT = 0.1;
    return (fabs(factor - 1.0) > LIMIT) ||
           (rAngle > LIMIT * 2*M_PI / TimeStep::stepsPerYear);
}


// Every TimeStep::interval days:
void FixedEmergence::update () {
    if (TimeStep::simulation > larvicidingEndStep) {
        larvicidingEndStep = TimeStep::future;
        larvicidingIneffectiveness = 1.0;
    }
}

// -----  Summary and intervention functions  -----

void FixedEmergence::intervLarviciding (const scnXml::LarvicidingDescAnoph& elt) {
    larvicidingIneffectiveness = 1 - elt.getEffectiveness().getValue();
    larvicidingEndStep = TimeStep::simulation + TimeStep(elt.getDuration().getValue());
}

double FixedEmergence::getLastN_v0 () const{
    double timestep_N_v0 = 0.0;
    int firstDay = TimeStep::simulation.inDays() - TimeStep::interval + 1;
    for (size_t i = 0; i < (size_t)TimeStep::interval; ++i) {
        size_t dYear1 = (firstDay + i - 1) % TimeStep::fromYears(1).inDays();
        timestep_N_v0 += mosqEmergeRate[dYear1];
    }
    return timestep_N_v0;
}

}
}
}
